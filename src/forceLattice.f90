    !********************************************************************
    !     DEMBody 3.0
    !     ***********
    !
    !     Force for all Parallel Lattice.
    !     -------------------------------
    !     @Using parallel lattice instead of particle list.
    !
    !     @Using rolling friction similar to LIGGGHTS
    !     @Using Hertz-Mindlin contact model similar to Wada
    !     @Using Cohesion model similar to Scheeres
    !     @Using damping model applicable for ice ball
    !
    !********************************************************************
    subroutine forceLattice()

    use global
    implicit none

    real(kind=8)  Dist(3),DistS,DistL,DistR,DistU(3)
    real(kind=8)  Vrel(3),Vrot(3),Vtot(3),ERR,Vnor(3),Vtan(3)
    real(kind=8)  normal_force(3),normal_forceL
    real(kind=8)  tangential_force(3),tangential_forceL
    real(kind=8)  rolling_moment(3),rolling_momentL
    real(kind=8)  cohesive_force(3)
    real(kind=8)  Ap,An
    real(kind=8)  Rij,Mij,Iij
    real(kind=8)  Kn,Cn,Kt,Ct,Kr,Cr,lnCOR
    real(kind=8)  Dn,Ds(3),DsL,Dtheta(3),DthetaL
    real(kind=8)  H(3),M(3)
    real(kind=8)  RV(3)
    logical :: slipping,rolling             !  State of friction of T
    logical :: touching                     !  State of touch of T
    integer :: I,J,K,L,LenNode              !  Linklist
    integer :: num1,num2,limit              !  Neighbor
    type(Nodelink),pointer :: Temp          !  Temporary pointer
    type(Nodelink),pointer :: TempH         !  Contact pointer
    real(kind=8) ShearPBC                   !  shearPBC
    real(kind=8) accMag                     !  Magtitude of contact acceleration                
    integer :: region                       !  Upper limit of Loops 
    logical :: check                        !  Check for meshgrid

    integer :: II                           !  Parallel Lattice
    integer :: NoInner,NoOuter              !  Number of Inner & Outer Particles
    integer :: IDInner(NLAT),IDOuter(NLAT)  !  ID array of Inner & Outer Particles
    integer :: particleI,particleJ          !  ID of two interacting particle

    F = 0.0D0
    FM = 0.0D0
    Energy = 0.0D0

    !  Loop over all bodies through the NodeTree.
    !$OMP PARALLEL DO REDUCTION(+:F) REDUCTION(+:FM) &
    !$OMP REDUCTION(+:Energy) &
    !$OMP& PRIVATE(check,shearPBC,Temp,TempH,LenNode,num1,num2,limit,&
    !$OMP& I,J,K,L,Dist,DistS,DistL,DistR,DistU,Vrel,Vrot,Vtot,ERR,Vnor,Vtan,&
    !$OMP& normal_force,normal_forceL,tangential_force,tangential_forceL,&
    !$OMP& rolling_moment,rolling_momentL,cohesive_force,Ap,An,Rij,Mij,Iij,&
    !$OMP& Kn,Cn,Kt,Ct,Kr,Cr,lnCOR,Dn,Ds,DsL,Dtheta,DthetaL,H,M,RV,&
    !$OMP& slipping,rolling,touching,&
    !$OMP& II,NoInner,NoOuter,IDInner,IDOuter,particleI,particleJ) SCHEDULE(DYNAMIC)
    !  Loop over all lattices
    do II = 1,LatNum
        NoInner = DEM(II)%NoInner
        NoOuter = DEM(II)%NoOuter
        IDInner = DEM(II)%IDInner
        IDOuter = DEM(II)%IDOuter

        if (NoInner .NE. 0) then
            !  Loop over inner region
            do I = 1,NoInner
                do J = 1,NoInner

                particleI = IDInner(I)
                particleJ = IDInner(J)

                !  calculated the NodeTree to detect the overlaps.
                limit = 0
                num1 = Linklist(particleJ) - Linklist(particleI)
                !!!!!!!!!!!!!!!!!!!
                select case(num1)
                case(CASE1,CASE2,CASE3,CASE4,CASE5,CASE6,CASE7,CASE8,CASE9,CASE10,CASE11,CASE12,CASE13)
                    limit = 1
                case (0)
                    if (J .GT. I) then
                        limit = 1
                    end if
                end select

                check = limit.EQ.1

                if (check) then   !  .OR. J.EQ.N
                    !  Initialize state params
                    do K = 1,3
                        Dist(K) = X(K,particleJ) -  X(K,particleI)
                        H(K) = 0.0
                        M(K) = 0.0
                    end do
                    touching = .false.
                    slipping = .false.
                    rolling = .false.
                    !  Distance vector
                    DistS = Dist(1)*Dist(1) + Dist(2)*Dist(2) + Dist(3)*Dist(3)
                    DistL = sqrt(DistS)
                    Dn = R(particleI) + R(particleJ) - DistL
                    !  lookup this contact(or not) in the Hertz list
                    LenNode = Head(particleI)%No
                    Temp => Head(particleI)
                    if (LenNode .NE. 0) then
                        Temp => Head(particleI)%next
                        do L = 1,LenNode
                            if (Temp%No .EQ. particleJ) then
                                do K = 1,3
                                    H(K) = Temp%Hertz(K)
                                    M(K) = Temp%Mrot(K)
                                end do
                                touching = Temp%is_touching
                                slipping = Temp%is_slipping
                                rolling = Temp%is_rolling
                                exit
                            else if (Temp%No.LT.particleJ .AND. associated(Temp%next)) then
                                Temp => Temp%next
                            else if (Temp%No .GT. particleJ) then
                                Temp => Temp%prev
                                exit
                            end if
                        end do
                    end if
                    !  When collision calculate the repulsive restoring spring force which is generated along the normal and tangential according to Hooke's law
                    if (Dn .GT. 0.0) then
                        DistR = 1.0/DistL
                        !  calculate the normal vector
                        do K=1,3
                            DistU(K) = Dist(K)*DistR
                        end do
                        Ap = (R(particleI)*R(particleI)-R(particleJ)*R(particleJ)+DistS)/2.0*DistR
                        An = DistL-Ap
#ifdef HertzMindlinVisco               
                        !  calculate material constant
                        Rij = R(particleI)*R(particleJ)/(R(particleI)+R(particleJ))
                        Mij = Body(particleI)*Body(particleJ)/(Body(particleI)+Body(particleJ))
                        Iij = 3.5*Inertia(particleI)*Inertia(particleJ)/(Inertia(particleI)+Inertia(particleJ))
                        Kn = 2.0*m_E*sqrt(Rij)/(3.0*(1.0-m_nu*m_nu))
                        Cn = -Kn*m_A*sqrt(Dn)
                        Kt = 2.0*m_E/(1.0+m_nu)/(2.0-m_nu)*sqrt(Rij)*sqrt(Dn)
                        !  select tangential damping mode
                        if (m_COR > 1.0) then
                            Ct = -2.0*m_E/(1.0+m_nu)/(2.0-m_nu)*sqrt(Dn)*m_A
                        elseif (m_COR >= 0.0) then
                            lnCOR=log(m_COR)
                            Ct = 2.0*sqrt(5.0/6.0)*lnCOR/sqrt(lnCOR**2+3.1415926**2) &
                            & *sqrt(2.0*Mij*m_E/(1.0+m_nu)/(2.0-m_nu))*(Rij**0.25)*(Dn**0.25)
                        else
                            Ct = 0
                        end if
                        Kr = 2.25*(Rij**2)*(m_mu_r**2)*Kn*sqrt(Dn)
                        Cr = 2.0*m_nita_r*sqrt(Iij*Kr)
#elif HertzMindlinResti
                        !  calculate material constant
                        Rij = R(particleI)*R(particleJ)/(R(particleI)+R(particleJ))
                        Mij = Body(particleI)*Body(particleJ)/(Body(particleI)+Body(particleJ))
                        Iij = 3.5*Inertia(particleI)*Inertia(particleJ)/(Inertia(particleI)+Inertia(particleJ))
                        Kn = 2.0*m_E*sqrt(Rij)/(3.0*(1.0-m_nu*m_nu))
                        lnCOR = log(m_COR)
                        Cn = 2.0*sqrt(5.0/6.0)*lnCOR/sqrt(lnCOR**2+3.1415926**2) &
                        & *sqrt(Mij*m_E/(1.0-m_nu*m_nu))*(Rij**0.25)*(Dn**0.25)
                        Kt = 2.0*m_E/(1.0+m_nu)/(2.0-m_nu)*sqrt(Rij)*sqrt(Dn)
                        Ct = 2.0*sqrt(5.0/6.0)*lnCOR/sqrt(lnCOR**2+3.1415926**2) &
                        & *sqrt(2.0*Mij*m_E/(1.0+m_nu)/(2.0-m_nu))*(Rij**0.25)*(Dn**0.25)
                        Kr = 2.25*(Rij**2)*(m_mu_r**2)*Kn*sqrt(Dn)
                        Cr = 2.0*m_nita_r*sqrt(Iij*Kr)
#endif
                        !  translate relative velocity
                        do K = 1,3
                            Vrel(K) = Xdot(K,particleJ) - Xdot(K,particleI)
                        end do
                        !  negative rotate relative velocity
                        Vrot(1) = (DistU(2)*W(3,particleJ) - DistU(3)*W(2,particleJ))*An + (DistU(2)*W(3,particleI) - DistU(3)*W(2,particleI))*Ap
                        Vrot(2) = (DistU(3)*W(1,particleJ) - DistU(1)*W(3,particleJ))*An + (DistU(3)*W(1,particleI) - DistU(1)*W(3,particleI))*Ap
                        Vrot(3) = (DistU(1)*W(2,particleJ) - DistU(2)*W(1,particleJ))*An + (DistU(1)*W(2,particleI) - DistU(2)*W(1,particleI))*Ap
                        !  totle relative velocity
                        do K = 1,3
                            Vtot(K) = Vrel(K) + Vrot(K)
                        end do
                        !  calculate the normal components of velocity
                        ERR = DistU(1)*Vtot(1) + DistU(2)*Vtot(2) + DistU(3)*Vtot(3)
                        do K = 1,3
                            Vnor(K) = ERR * DistU(K)
                        end do
                        !  tangential velocity
                        do K = 1,3
                            Vtan(K) = Vtot(K) - Vnor(K)
                        end do

                        !  normal force of Particle J
                        do K = 1,3
                            normal_force(K) = Kn*(Dn**1.5)*DistU(K) + Cn*Vnor(K)
                        end do
                        normal_forceL = sqrt(normal_force(1)*normal_force(1) + normal_force(2)*normal_force(2) + normal_force(3)*normal_force(3))

                        !  Add energy
                        Energy(particleI) = Energy(particleI) + 0.2*Kn*(Dn**2.5)
                        Energy(particleJ) = Energy(particleJ) + 0.2*Kn*(Dn**2.5)

                        !  tangential deform
                        do K = 1,3
                            Ds(K) = Vtan(K)*Dt
                        end do
                        DsL = sqrt(Ds(1)*Ds(1) + Ds(2)*Ds(2) + Ds(3)*Ds(3))

                        !  tangential force of Particle J
                        do K = 1,3
                            tangential_force(K) = - Kt*Ds(K) + Ct*Vtan(K) + H(K)
                        end do
                        tangential_forceL = sqrt(tangential_force(1)*tangential_force(1) + tangential_force(2)*tangential_force(2) + tangential_force(3)*tangential_force(3))

                        if (slipping) then  !  Have slipped
                            if (DsL .GT. 1e-8) then  !  Still slipping
                                do K = 1,3
                                    tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    tangential_force(K) = 0.0  !  Particle J
                                end do
                                slipping = .false.
                            end if
                        else
                            if (tangential_forceL .GT. normal_forceL*m_mu_s) then  !  Slipping
                                slipping = .true.
                                if (DsL .GT. 1e-14) then
                                    do K = 1,3
                                        tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL
                                    end do
                                else
                                    do K = 1,3
                                        tangential_force(K) = 0.0
                                    end do
                                end if
                            else
                                slipping = .false.  !  Sticking
                            end if
                        end if

                        touching = .true.
                        !  Apply force
                        do K = 1,3
                            F(K,particleJ) = + normal_force(K) + tangential_force(K) + F(K,particleJ)
                            F(K,particleI) = - normal_force(K) - tangential_force(K) + F(K,particleI)
                        end do
                        !  Apply moment
                        FM(1,particleJ) = - An*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,particleJ) 
                        FM(2,particleJ) = - An*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,particleJ) 
                        FM(3,particleJ) = - An*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,particleJ)                     
                        FM(1,particleI) = - Ap*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,particleI) 
                        FM(2,particleI) = - Ap*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,particleI) 
                        FM(3,particleI) = - Ap*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,particleI)

                        !  rolling deform
                        do K = 1,3
                            Dtheta(K) = (W(K,particleI)-W(K,particleJ))*Dt
                        end do
                        DthetaL = sqrt(Dtheta(1)*Dtheta(1) + Dtheta(2)*Dtheta(2) + Dtheta(3)*Dtheta(3))

                        !  rolling moment of Particle J
                        do K = 1,3
                            rolling_moment(K) = Kr*Dtheta(K) + M(K)
                        end do
                        rolling_momentL = sqrt(rolling_moment(1)*rolling_moment(1) + rolling_moment(2)*rolling_moment(2) + rolling_moment(3)*rolling_moment(3))                

                        if (rolling) then  !  Have rolled
                            if (DthetaL .GT. 1e-8) then  !  Still slipping
                                do K = 1,3
                                    rolling_moment(K) = m_mu_r*Rij*normal_forceL*Dtheta(K)/DthetaL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    rolling_moment(K) = 0.0  !  Particle J
                                end do
                                rolling = .false.
                            end if
                        else
                            if (rolling_momentL .GT. normal_forceL*Rij*m_mu_r) then  !  Rolling
                                rolling = .true.
                                if (DthetaL .GT. 1e-14) then
                                    do K = 1,3
                                        rolling_moment(K) = m_mu_r*Rij*normal_forceL*Dtheta(K)/DthetaL
                                    end do
                                else
                                    do K = 1,3
                                        rolling_moment(K) = 0.0
                                    end do
                                end if
                            else
                                rolling = .false.  !  Sticking
                                do K = 1,3
                                    rolling_moment(K) = rolling_moment(K) + Cr*(W(K,particleI)-W(K,particleJ))
                                end do
                            end if
                        end if                                  
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
                        !  Apply moment
                        do K = 1,3
                            FM(K,particleJ) = + rolling_moment(K) + FM(K,particleJ) 
                            FM(K,particleI) = - rolling_moment(K) + FM(K,particleI) 
                        end do                                          
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

                        !  cohesive force
                        do K = 1,3
                            cohesive_force(K) = - 1.0*m_Beta*3.1415926*Rij**2*DistU(K)
                        end do
                        do K = 1,3
                            F(K,particleJ) = F(K,particleJ) + cohesive_force(K)
                            F(K,particleI) = F(K,particleI) - cohesive_force(K)
                        end do
                        !  memory the contact in the Hertz linklist.
                        if (associated(Temp%prev)) then
                            !  Temp is in center of linklist!!!
                            if (Temp%No .EQ. particleJ) then
                                !  Have contacted.
                                do K = 1,3
                                    Temp%Hertz(K) = tangential_force(K)
                                    Temp%Mrot(K) = rolling_moment(K)
                                end do
                                Temp%is_touching = touching
                                Temp%is_slipping = slipping
                                Temp%is_rolling = rolling
                            else
                                !  First contacted.
                                allocate(TempH)
                                TempH = Nodelink(particleJ,tangential_force,rolling_moment,&
                                & touching,slipping,rolling,Temp,Temp%next)
                                if (associated(Temp%next)) Temp%next%prev => TempH
                                Temp%next => TempH
                                Head(particleI)%No = LenNode + 1
                            end if
                        else
                            !  Temp is Head of linklist!!!
                            allocate(TempH)
                            TempH = Nodelink(particleJ,tangential_force,rolling_moment,&
                            & touching,slipping,rolling,Temp,Temp%next)
                            if (associated(Temp%next)) Temp%next%prev => TempH
                            Temp%next => TempH
                            Head(particleI)%No = LenNode + 1
                        end if
                    else
                        !  memory the separation in the Hertz linklist.
                        if (associated(Temp%prev) .AND. Temp%No.EQ.particleJ) then
                            !  Temp is center of linklist!!!
                            Temp%prev%next => Temp%next
                            if(associated(Temp%next)) Temp%next%prev => Temp%prev
                            Head(particleI)%No = LenNode - 1
                            deallocate(Temp)
                            !  When else Temp is Head of linklist!!!
                        end if

                        Rij = R(particleI)*R(particleJ)/(R(particleI)+R(particleJ))
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        if (Dn > -Rij*(m_r_cut-1.0)*0.5) then
                            DistR = 1.0/DistL
                            !  calculate the normal vector
                            do K=1,3
                                DistU(K) = Dist(K)*DistR
                            end do
                            !  cohesive force
                            do K = 1,3
                                cohesive_force(K) = - 1.0*m_Beta*3.1415926*Rij**2*DistU(K)
                            end do
                            do K = 1,3
                                F(K,particleJ) = F(K,particleJ) + cohesive_force(K)
                                F(K,particleI) = F(K,particleI) - cohesive_force(K)
                            end do                      
                        else if (Dn > -Rij*(m_r_cut-1.0)) then
                            DistR = 1.0/DistL
                            !  calculate the normal vector
                            do K=1,3
                                DistU(K) = Dist(K)*DistR
                            end do
                            !  cohesive force
                            do K = 1,3
                                cohesive_force(K) = - 1.0*m_Beta*3.1415926*Rij**2 &
                                &*2.0*(Dn/(R(particleI)+R(particleJ))/(m_r_cut-1) + 1.0)*DistU(K)
                            end do
                            do K = 1,3
                                F(K,particleJ) = F(K,particleJ) + cohesive_force(K)
                                F(K,particleI) = F(K,particleI) - cohesive_force(K)
                            end do                                                   
                        end if
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    end if
                end if
                end do    
            end do
        end if

        if (NoInner.NE.0 .AND. NoOuter.NE.0) then
            !  Loop over outer region
            do I = 1,NoInner
                do J = 1,NoOuter

                particleI = IDInner(I)
                particleJ = IDOuter(J)

                !  calculated the NodeTree to detect the overlaps.
                limit = 0
                num1 = ABS(Linklist(particleJ) - Linklist(particleI))
                !!!!!!!!!!!!!!!!!!!
                select case(num1)
                case(CASE1,-CASE2,-CASE3,-CASE4,-CASE5,-CASE6,-CASE7,-CASE8,-CASE9,-CASE10,-CASE11,-CASE12,-CASE13)
                    limit = 1
                case (0)
                    limit = 1
                end select

                check = limit.EQ.1

                if (check) then   !  .OR. J.EQ.N
                    !  Initialize state params
                    do K = 1,3
                        Dist(K) = X(K,particleJ) -  X(K,particleI)
                        H(K) = 0.0
                        M(K) = 0.0
                    end do
                    touching = .false.
                    slipping = .false.
                    rolling = .false.
                    !  Distance vector
                    DistS = Dist(1)*Dist(1) + Dist(2)*Dist(2) + Dist(3)*Dist(3)
                    DistL = sqrt(DistS)
                    Dn = R(particleI) + R(particleJ) - DistL
                    !  lookup this contact(or not) in the Hertz list
                    LenNode = Head(particleI)%No
                    Temp => Head(particleI)
                    if (LenNode .NE. 0) then
                        Temp => Head(particleI)%next
                        do L = 1,LenNode
                            if (Temp%No .EQ. particleJ) then
                                do K = 1,3
                                    H(K) = Temp%Hertz(K)
                                    M(K) = Temp%Mrot(K)
                                end do
                                touching = Temp%is_touching
                                slipping = Temp%is_slipping
                                rolling = Temp%is_rolling
                                exit
                            else if (Temp%No.LT.particleJ .AND. associated(Temp%next)) then
                                Temp => Temp%next
                            else if (Temp%No .GT. particleJ) then
                                Temp => Temp%prev
                                exit
                            end if
                        end do
                    end if
                    !  When collision calculate the repulsive restoring spring force which is generated along the normal and tangential according to Hooke's law
                    if (Dn .GT. 0.0) then
                        DistR = 1.0/DistL
                        !  calculate the normal vector
                        do K=1,3
                            DistU(K) = Dist(K)*DistR
                        end do
                        Ap = (R(particleI)*R(particleI)-R(particleJ)*R(particleJ)+DistS)/2.0*DistR
                        An = DistL-Ap
#ifdef HertzMindlinVisco                    
                        !  calculate material constant
                        Rij = R(particleI)*R(particleJ)/(R(particleI)+R(particleJ))
                        Mij = Body(particleI)*Body(particleJ)/(Body(particleI)+Body(particleJ))
                        Iij = 3.5*Inertia(particleI)*Inertia(particleJ)/(Inertia(particleI)+Inertia(particleJ))
                        Kn = 2.0*m_E*sqrt(Rij)/(3.0*(1.0-m_nu*m_nu))
                        Cn = -Kn*m_A*sqrt(Dn)
                        Kt = 2.0*m_E/(1.0+m_nu)/(2.0-m_nu)*sqrt(Rij)*sqrt(Dn)
                        !  select tangential damping mode
                        if (m_COR > 1.0) then
                            Ct = -2.0*m_E/(1.0+m_nu)/(2.0-m_nu)*sqrt(Dn)*m_A
                        elseif (m_COR >= 0.0) then
                            lnCOR=log(m_COR)
                            Ct = 2.0*sqrt(5.0/6.0)*lnCOR/sqrt(lnCOR**2+3.1415926**2) &
                            & *sqrt(2.0*Mij*m_E/(1.0+m_nu)/(2.0-m_nu))*(Rij**0.25)*(Dn**0.25)
                        else
                            Ct = 0
                        end if
                        Kr = 2.25*(Rij**2)*(m_mu_r**2)*Kn*sqrt(Dn)
                        Cr = 2.0*m_nita_r*sqrt(Iij*Kr)
#elif HertzMindlinResti
                        !  calculate material constant
                        Rij = R(particleI)*R(particleJ)/(R(particleI)+R(particleJ))
                        Mij = Body(particleI)*Body(particleJ)/(Body(particleI)+Body(particleJ))
                        Iij = 3.5*Inertia(particleI)*Inertia(particleJ)/(Inertia(particleI)+Inertia(particleJ))
                        Kn = 2.0*m_E*sqrt(Rij)/(3.0*(1.0-m_nu*m_nu))
                        lnCOR = log(m_COR)
                        Cn = 2.0*sqrt(5.0/6.0)*lnCOR/sqrt(lnCOR**2+3.1415926**2) &
                        & *sqrt(Mij*m_E/(1.0-m_nu*m_nu))*(Rij**0.25)*(Dn**0.25)
                        Kt = 2.0*m_E/(1.0+m_nu)/(2.0-m_nu)*sqrt(Rij)*sqrt(Dn)
                        Ct = 2.0*sqrt(5.0/6.0)*lnCOR/sqrt(lnCOR**2+3.1415926**2) &
                        & *sqrt(2.0*Mij*m_E/(1.0+m_nu)/(2.0-m_nu))*(Rij**0.25)*(Dn**0.25)
                        Kr = 2.25*(Rij**2)*(m_mu_r**2)*Kn*sqrt(Dn)
                        Cr = 2.0*m_nita_r*sqrt(Iij*Kr)
#endif
                        !  translate relative velocity
                        do K = 1,3
                            Vrel(K) = Xdot(K,particleJ) - Xdot(K,particleI)
                        end do
                        !  negative rotate relative velocity
                        Vrot(1) = (DistU(2)*W(3,particleJ) - DistU(3)*W(2,particleJ))*An + (DistU(2)*W(3,particleI) - DistU(3)*W(2,particleI))*Ap
                        Vrot(2) = (DistU(3)*W(1,particleJ) - DistU(1)*W(3,particleJ))*An + (DistU(3)*W(1,particleI) - DistU(1)*W(3,particleI))*Ap
                        Vrot(3) = (DistU(1)*W(2,particleJ) - DistU(2)*W(1,particleJ))*An + (DistU(1)*W(2,particleI) - DistU(2)*W(1,particleI))*Ap
                        !  totle relative velocity
                        do K = 1,3
                            Vtot(K) = Vrel(K) + Vrot(K)
                        end do
                        !  calculate the normal components of velocity
                        ERR = DistU(1)*Vtot(1) + DistU(2)*Vtot(2) + DistU(3)*Vtot(3)
                        do K = 1,3
                            Vnor(K) = ERR * DistU(K)
                        end do
                        !  tangential velocity
                        do K = 1,3
                            Vtan(K) = Vtot(K) - Vnor(K)
                        end do

                        !  normal force of Particle J
                        do K = 1,3
                            normal_force(K) = Kn*(Dn**1.5)*DistU(K) + Cn*Vnor(K)
                        end do
                        normal_forceL = sqrt(normal_force(1)*normal_force(1) + normal_force(2)*normal_force(2) + normal_force(3)*normal_force(3))

                        !  Add energy
                        Energy(particleI) = Energy(particleI) + 0.2*Kn*(Dn**2.5)
                        Energy(particleJ) = Energy(particleJ) + 0.2*Kn*(Dn**2.5)

                        !  tangential deform
                        do K = 1,3
                            Ds(K) = Vtan(K)*Dt
                        end do
                        DsL = sqrt(Ds(1)*Ds(1) + Ds(2)*Ds(2) + Ds(3)*Ds(3))

                        !  tangential force of Particle J
                        do K = 1,3
                            tangential_force(K) = - Kt*Ds(K) + Ct*Vtan(K) + H(K)
                        end do
                        tangential_forceL = sqrt(tangential_force(1)*tangential_force(1) + tangential_force(2)*tangential_force(2) + tangential_force(3)*tangential_force(3))

                        if (slipping) then  !  Have slipped
                            if (DsL .GT. 1e-8) then  !  Still slipping
                                do K = 1,3
                                    tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    tangential_force(K) = 0.0  !  Particle J
                                end do
                                slipping = .false.
                            end if
                        else
                            if (tangential_forceL .GT. normal_forceL*m_mu_s) then  !  Slipping
                                slipping = .true.
                                if (DsL .GT. 1e-14) then
                                    do K = 1,3
                                        tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL
                                    end do
                                else
                                    do K = 1,3
                                        tangential_force(K) = 0.0
                                    end do
                                end if
                            else
                                slipping = .false.  !  Sticking
                            end if
                        end if

                        touching = .true.
                        !  Apply force
                        do K = 1,3
                            F(K,particleJ) = + normal_force(K) + tangential_force(K) + F(K,particleJ)
                            F(K,particleI) = - normal_force(K) - tangential_force(K) + F(K,particleI)
                        end do
                        !  Apply moment
                        FM(1,particleJ) = - An*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,particleJ) 
                        FM(2,particleJ) = - An*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,particleJ) 
                        FM(3,particleJ) = - An*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,particleJ)                     
                        FM(1,particleI) = - Ap*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,particleI) 
                        FM(2,particleI) = - Ap*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,particleI) 
                        FM(3,particleI) = - Ap*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,particleI)

                        !  rolling deform
                        do K = 1,3
                            Dtheta(K) = (W(K,particleI)-W(K,particleJ))*Dt
                        end do
                        DthetaL = sqrt(Dtheta(1)*Dtheta(1) + Dtheta(2)*Dtheta(2) + Dtheta(3)*Dtheta(3))

                        !  rolling moment of Particle J
                        do K = 1,3
                            rolling_moment(K) = Kr*Dtheta(K) + M(K)
                        end do
                        rolling_momentL = sqrt(rolling_moment(1)*rolling_moment(1) + rolling_moment(2)*rolling_moment(2) + rolling_moment(3)*rolling_moment(3))                

                        if (rolling) then  !  Have rolled
                            if (DthetaL .GT. 1e-8) then  !  Still slipping
                                do K = 1,3
                                    rolling_moment(K) = m_mu_r*Rij*normal_forceL*Dtheta(K)/DthetaL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    rolling_moment(K) = 0.0  !  Particle J
                                end do
                                rolling = .false.
                            end if
                        else
                            if (rolling_momentL .GT. normal_forceL*Rij*m_mu_r) then  !  Rolling
                                rolling = .true.
                                if (DthetaL .GT. 1e-14) then
                                    do K = 1,3
                                        rolling_moment(K) = m_mu_r*Rij*normal_forceL*Dtheta(K)/DthetaL
                                    end do
                                else
                                    do K = 1,3
                                        rolling_moment(K) = 0.0
                                    end do
                                end if
                            else
                                rolling = .false.  !  Sticking
                                do K = 1,3
                                    rolling_moment(K) = rolling_moment(K) + Cr*(W(K,particleI)-W(K,particleJ))
                                end do
                            end if
                        end if                                  
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
                        !  Apply moment
                        do K = 1,3
                            FM(K,particleJ) = + rolling_moment(K) + FM(K,particleJ) 
                            FM(K,particleI) = - rolling_moment(K) + FM(K,particleI) 
                        end do                                          
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

                        !  cohesive force
                        do K = 1,3
                            cohesive_force(K) = - 1.0*m_Beta*3.1415926*Rij**2*DistU(K)
                        end do
                        do K = 1,3
                            F(K,particleJ) = F(K,particleJ) + cohesive_force(K)
                            F(K,particleI) = F(K,particleI) - cohesive_force(K)
                        end do
                        !  memory the contact in the Hertz linklist.
                        if (associated(Temp%prev)) then
                            !  Temp is in center of linklist!!!
                            if (Temp%No .EQ. particleJ) then
                                !  Have contacted.
                                do K = 1,3
                                    Temp%Hertz(K) = tangential_force(K)
                                    Temp%Mrot(K) = rolling_moment(K)
                                end do
                                Temp%is_touching = touching
                                Temp%is_slipping = slipping
                                Temp%is_rolling = rolling
                            else
                                !  First contacted.
                                allocate(TempH)
                                TempH = Nodelink(particleJ,tangential_force,rolling_moment,&
                                & touching,slipping,rolling,Temp,Temp%next)
                                if (associated(Temp%next)) Temp%next%prev => TempH
                                Temp%next => TempH
                                Head(particleI)%No = LenNode + 1
                            end if
                        else
                            !  Temp is Head of linklist!!!
                            allocate(TempH)
                            TempH = Nodelink(particleJ,tangential_force,rolling_moment,&
                            & touching,slipping,rolling,Temp,Temp%next)
                            if (associated(Temp%next)) Temp%next%prev => TempH
                            Temp%next => TempH
                            Head(particleI)%No = LenNode + 1
                        end if
                    else
                        !  memory the separation in the Hertz linklist.
                        if (associated(Temp%prev) .AND. Temp%No.EQ.particleJ) then
                            !  Temp is center of linklist!!!
                            Temp%prev%next => Temp%next
                            if(associated(Temp%next)) Temp%next%prev => Temp%prev
                            Head(particleI)%No = LenNode - 1
                            deallocate(Temp)
                            !  When else Temp is Head of linklist!!!
                        end if

                        Rij = R(particleI)*R(particleJ)/(R(particleI)+R(particleJ))
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        if (Dn > -Rij*(m_r_cut-1.0)*0.5) then
                            DistR = 1.0/DistL
                            !  calculate the normal vector
                            do K=1,3
                                DistU(K) = Dist(K)*DistR
                            end do
                            !  cohesive force
                            do K = 1,3
                                cohesive_force(K) = - 1.0*m_Beta*3.1415926*Rij**2*DistU(K)
                            end do
                            do K = 1,3
                                F(K,particleJ) = F(K,particleJ) + cohesive_force(K)
                                F(K,particleI) = F(K,particleI) - cohesive_force(K)
                            end do                      
                        else if (Dn > -Rij*(m_r_cut-1.0)) then
                            DistR = 1.0/DistL
                            !  calculate the normal vector
                            do K=1,3
                                DistU(K) = Dist(K)*DistR
                            end do
                            !  cohesive force
                            do K = 1,3
                                cohesive_force(K) = - 1.0*m_Beta*3.1415926*Rij**2 &
                                &*2.0*(Dn/(R(particleI)+R(particleJ))/(m_r_cut-1) + 1.0)*DistU(K)
                            end do
                            do K = 1,3
                                F(K,particleJ) = F(K,particleJ) + cohesive_force(K)
                                F(K,particleI) = F(K,particleI) - cohesive_force(K)
                            end do                                                   
                        end if
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    end if
                end if
                end do    
            end do
        end if
    end do
    !$OMP END PARALLEL DO

    !  calculate mirrored force if using PBC & Shear PBC
    if (isPeriodic) then
        call forceMirror
    end if
    
    !  calculate force of contactable walls if using contactable walls
    if (isContactWall) then
        call forceContactWalls
    end if
    
    !  calculate force of bonded walls if using bonded walls
    if (isBondedWall) then                     
        call forceBondedWalls 
    end if
    
    !  calculate force of funnel walls if using funnel walls
    if (isFunnelWall) then
        call forceFunnelWalls
    end if

    !$OMP PARALLEL DO PRIVATE(I,K,accMag)
    do I = 1,N
        do K = 1,3
            F(K,I) = F(K,I)/Body(I)
            FM(K,I) = FM(K,I)/Inertia(I)
        end do
        accMag = sqrt(F(1,I)**2 + F(2,I)**2+F(3,I)**2)
        if (accMag .GE. MAX_ACC) then
            write(*,*) I,accMag
            do K = 1,3
                F(K,I) = F(K,I)/accMag*MAX_ACC
            end do
        end if
        F(3,I) = F(3,I) - G
    end do
    !$OMP END PARALLEL DO

    !  calculate Planet force if in Planet system
    if (isPlanet) then
        call Planet
    end if

    return
    end
