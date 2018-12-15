    !********************************************************************
    !     DEMBody 5.0
    !     ***********
    !
    !     Force for all Particles but using Parallel Lattice.
    !     -------------------------------
    !     @Using particle list coupled with parallel lattice.
    !     @This loop manner may be suit for the case in which the number
    !     @of particles is less than the number of lattices, or the case in
    !     @which the particles distribute nonuniformly (which means that many
    !     @lattices has no particle.) 
    !
    !     @Using rolling friction similar to LIGGGHTS
    !     @Using Hertz-Mindlin contact model similar to Wada
    !     @Using Cohesion model similar to Scheeres
    !     @Using damping model applicable for ice ball
    !
    !********************************************************************
    subroutine forceParticleLattice()

    use global
    use omp_lib
    implicit none

    real(8) ::  Dist(3),DistS,DistL,DistR,DistU(3)
    real(8) ::  Vrel(3),Vrot(3),Vtot(3),ERR,Vnor(3),Vtan(3)
    real(8) ::  normal_force(3),normal_forceL
    real(8) ::  tangential_force(3),tangential_forceL,tangential_history(3)
    real(8) ::  rolling_moment(3),rolling_momentL,rolling_history(3)
    real(8) ::  twisting_moment(3),twisting_momentL,twisting_history(3)
    real(8) ::  cohesive_force(3)
    real(8) ::  gravity_force(3)
    real(8) ::  Ap,An
    real(8) ::  Rij,Mij,Iij
    real(8) ::  Kn,Cn,Ks,Cs,Kr,Cr,Kt,Ct,lnCOR
    real(8) ::  Dn,Ds(3),DsL,Dtheta(3),DthetaL,DthetaR(3),DthetaRL,DthetaT(3),DthetaTL
    real(8) ::  H(3),Mr(3),Mt(3)
    real(8) ::  RV(3)
    logical :: slipping,rolling,twisting    !  State of friction of T
    logical :: touching                     !  State of touch of T
    integer :: I,J,K,L,LenNodeJ,LenNode     !  Linklist
    integer :: num1,num2,limit              !  Neighbor
    type(Nodelink),pointer :: Temp          !  Temporary pointer
    type(Nodelink),pointer :: TempH         !  Contact pointer
    type(Nodelink),pointer :: Tail          !  Tail pointer
    real(8) ShearPBC                        !  shearPBC
    real(8) accMag                          !  Magtitude of contact acceleration                
    integer :: region                       !  Upper limit of Loops 
    logical :: check                        !  Check for meshgrid

    integer :: JJ                           !  Parallel Lattice
    integer :: IDJ                          !  Parallel Lattice J (outer lattice)
    type(Neighbor),pointer :: IDInnerI      !  List of Inner Particles
    type(Neighbor),pointer :: IDInnerJ      !  List of Inner Particles
    type(Neighbor),pointer :: IDOuterJ      !  List of Outer Particles   
    integer :: particleI,particleJ          !  ID of two interacting particle
    
    real(8) :: ostart,oend

    character(30) :: FileNameHead
    character(30) :: FileNameForce
    
    F = 0.0D0
    FM = 0.0D0
    Energy = 0.0D0
    Heat = 0.0D0
    
    !write(FileNameHead,'(F10.5)') Time
    !FileNameHead = trim(FileNameHead)//'Head.txt'
    !open(123,FILE=FileNameHead)
    
    !ostart = omp_get_wtime()
    !  Loop over all bodies through the NodeTree.
    !$OMP PARALLEL DO &
    !$OMP& PRIVATE(check,shearPBC,LenNodeJ,LenNode,Temp,TempH,Tail,num1,num2,limit,&
    !$OMP& I,J,K,L,Dist,DistS,DistL,DistR,DistU,Vrel,Vrot,Vtot,ERR,Vnor,Vtan,&
    !$OMP& normal_force,normal_forceL,tangential_force,tangential_forceL,tangential_history,&
    !$OMP& rolling_moment,rolling_momentL,rolling_history,twisting_moment,twisting_momentL,twisting_history,cohesive_force,gravity_force,Ap,An,Rij,Mij,Iij,&
    !$OMP& Kn,Cn,Ks,Cs,Kr,Cr,Kt,Ct,lnCOR,Dn,Ds,DsL,Dtheta,DthetaL,DthetaR,DthetaRL,DthetaT,DthetaTL,H,Mr,Mt,RV,&
    !$OMP& slipping,rolling,twisting,touching,&
    !$OMP& JJ,IDJ,IDInnerI,IDInnerJ,IDOuterJ,particleI,particleJ) SCHEDULE(GUIDED)
    ! Loop over all lattices
    do I = 1,N
        
        particleI = I
        Tail => Head(particleI)

        !  Loop over inner region
        IDInnerJ => IDInner(Linklist(I))
        do while (associated(IDInnerJ%next))
            IDInnerJ => IDInnerJ%next
            particleJ = IDInnerJ%No    

            !if (particleI.EQ.PP .OR. particleJ.EQ.PP) then
            !    write(123,'(F15.5,2X,A5,2X,2I5,2X)',advance='no') Time,'Inner',particleI,particleJ
            !end if

            !if (particleI.EQ.PP .OR. particleJ.EQ.PP) then
            !    write(123,'(A6,2X,I5,2X)',advance='no') 'check:',check
            !end if

            if (particleI .NE. particleJ) then
                !  Initialize state params
                do K = 1,3
                    Dist(K) = X(K,particleJ) -  X(K,particleI)
                    H(K) = 0.0D0
                    Mr(K) = 0.0D0
                    Mt(K) = 0.0D0
                end do
                touching = .false.
                slipping = .false.
                rolling = .false.
                twisting = .false.
                !  Distance vector
                DistS = Dist(1)*Dist(1) + Dist(2)*Dist(2) + Dist(3)*Dist(3)
                DistL = sqrt(DistS)
                Dn = R(particleI) + R(particleJ) - DistL
                !  lookup this contact(or not) in the Hertz list I
                LenNode = Head(particleI)%No
                Temp => Tail
                do while(associated(Temp%next))
                    Temp => Temp%next
                    if (Temp%No .EQ. particleJ) then
                        do K = 1,3
                            H(K) = Temp%Hertz(K)
                            Mr(K) = Temp%Mrot(K)
                            Mt(K) = Temp%Mtwist(K)
                        end do
                        touching = Temp%is_touching
                        slipping = Temp%is_slipping
                        rolling = Temp%is_rolling
                        twisting = Temp%is_twisting
                        exit
                    else if (Temp%No .GT. particleJ) then
                        Temp => Temp%prev
                        exit
                    end if
                end do      
                !if (particleI.EQ.PP .OR. particleJ.EQ.PP) then
                !    write(123,'(F15.5,2X)',advance='no') Dn
                !end if
                !  When collision calculate the repulsive restoring spring force which is generated along the normal and tangential according to Hooke's law
                if (Dn .GT. 0.0D0) then
                    DistR = 1.0D0/DistL
                    !  calculate the normal vector
                    do K=1,3
                        DistU(K) = Dist(K)*DistR
                    end do
                    Ap = (R(particleI)*R(particleI)-R(particleJ)*R(particleJ)+DistS)/2.0D0*DistR
                    An = DistL-Ap
#ifdef HertzMindlinVisco    
                    !  calculate material constant
                    Rij = R(particleI)*R(particleJ)/(R(particleI)+R(particleJ))
                    Mij = Body(particleI)*Body(particleJ)/(Body(particleI)+Body(particleJ))
                    Kn = 2.0D0*m_E*sqrt(Rij*Dn)/(3.0D0*(1.0D0-m_nu*m_nu))
                    Cn = -Kn*m_A
                    Ks = 2.0D0*m_E/(1.0D0+m_nu)/(2.0D0-m_nu)*sqrt(Rij)*sqrt(Dn)
                    !  select tangential damping mode
                    if (m_COR > 1.0D0) then
                        Cs = -2.0D0*m_E/(1.0D0+m_nu)/(2.0D0-m_nu)*sqrt(Dn)*m_A
                    elseif (m_COR >= 0.0D0) then
                        lnCOR=log(m_COR)
                        Cs = 2.0D0*sqrt(5.0D0/6.0D0)*lnCOR/sqrt(lnCOR**2+3.1415926D0**2) &
                        & *sqrt(2.0D0*Mij*m_E/(1.0D0+m_nu)/(2.0D0-m_nu))*(Rij**0.25)*(Dn**0.25)
                    else
                        Cs = 0.0D0
                    end if
                    Kr = 0.25D0*Kn*(m_Beta*Rij)**2
                    Cr = 0.25D0*Cn*(m_Beta*Rij)**2
                    Kt = 0.5D0*Ks*(m_Beta*Rij)**2
                    Ct = 0.5D0*Cs*(m_Beta*Rij)**2
#elif HertzMindlinResti
                    !  calculate material constant
                    Rij = R(particleI)*R(particleJ)/(R(particleI)+R(particleJ))
                    Mij = Body(particleI)*Body(particleJ)/(Body(particleI)+Body(particleJ))
                    Kn = 2.0D0*m_E*sqrt(Rij*Dn)/(3.0D0*(1.0D0-m_nu*m_nu))
                    lnCOR = log(m_COR)
                    Cn = 2.0D0*sqrt(5.0D0/6.0D0)*lnCOR/sqrt(lnCOR**2+3.1415926D0**2) &
                    & *sqrt(Mij*m_E/(1.0D0-m_nu*m_nu))*(Rij**0.25)*(Dn**0.25)
                    Ks = 2.0D0*m_E/(1.0D0+m_nu)/(2.0D0-m_nu)*sqrt(Rij)*sqrt(Dn)
                    Cs = 2.0D0*sqrt(5.0D0/6.0D0)*lnCOR/sqrt(lnCOR**2+3.1415926D0**2) &
                    & *sqrt(2.0D0*Mij*m_E/(1.0D0+m_nu)/(2.0D0-m_nu))*(Rij**0.25)*(Dn**0.25)
                    Kr = 0.25D0*Kn*(m_Beta*Rij)**2
                    Cr = 0.25D0*Cn*(m_Beta*Rij)**2
                    Kt = 0.5D0*Ks*(m_Beta*Rij)**2
                    Ct = 0.5D0*Cs*(m_Beta*Rij)**2
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
                        normal_force(K) = Kn*Dn*DistU(K) + Cn*Vnor(K)
                    end do
                    normal_forceL = sqrt(normal_force(1)*normal_force(1) + normal_force(2)*normal_force(2) + normal_force(3)*normal_force(3))

                    !  Add energy
                    Energy(particleI) = Energy(particleI) + 0.2D0*Kn*(Dn**2)
                    !Energy(particleJ) = Energy(particleJ) + 0.2D0*Kn*(Dn**2)

                    !  tangential deform
                    do K = 1,3
                        Ds(K) = Vtan(K)*Dt
                    end do
                    DsL = sqrt(Ds(1)*Ds(1) + Ds(2)*Ds(2) + Ds(3)*Ds(3))

                    !  tangential force of Particle J
                    do K = 1,3
                        tangential_force(K) = - Ks*Ds(K) + Cs*Vtan(K) + H(K)
                    end do
                    tangential_forceL = sqrt(tangential_force(1)*tangential_force(1) + tangential_force(2)*tangential_force(2) + tangential_force(3)*tangential_force(3))

                    !if (slipping) then  !  Have slipped
                    !    if (DsL .GT. 1.0e-14) then  !  Still slipping
                    !        do K = 1,3
                    !            tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL  !  Particle J
                    !        end do
                    !    else  !  Approach sticking
                    !        do K = 1,3
                    !            tangential_force(K) = 0.0D0  !  Particle J
                    !        end do
                    !        slipping = .false.
                    !    end if
                    !else
                        if (tangential_forceL .GT. normal_forceL*m_mu_s) then  !  Slipping
                            slipping = .true.
                            if (DsL .GT. 1.0e-14) then
                                do K = 1,3
                                    tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL
                                    tangential_history(K) = tangential_force(K)
                                end do
                                !  Add frictional heat
                                Heat(1,particleI) = Heat(1,particleI) + 0.5D0*m_mu_d*normal_forceL*DsL
                                !Heat(1,particleJ) = Heat(1,particleJ) + 0.5D0*m_mu_d*normal_forceL*DsL
                            else
                                do K = 1,3
                                    tangential_force(K) = 0.0D0
                                    tangential_history(K) = tangential_force(K)
                                end do
                            end if
                        else
                            slipping = .false.  !  Sticking
                            do K = 1,3
                                tangential_history(K) = tangential_force(K) - Cs*Vtan(K)
                            end do
                        end if
                    !end if

                    touching = .true.
                    !  Apply force
                    do K = 1,3
                        !F(K,particleJ) = + normal_force(K) + tangential_force(K) + F(K,particleJ)
                        F(K,particleI) = - normal_force(K) - tangential_force(K) + F(K,particleI)
                    end do
                    !  Apply moment
                    !FM(1,particleJ) = - An*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,particleJ) 
                    !FM(2,particleJ) = - An*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,particleJ) 
                    !FM(3,particleJ) = - An*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,particleJ)                     
                    FM(1,particleI) = - Ap*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,particleI) 
                    FM(2,particleI) = - Ap*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,particleI) 
                    FM(3,particleI) = - Ap*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,particleI)

                    !  rolling
                    do K = 1,3
                        Dtheta(K) = (W(K,particleI)-W(K,particleJ))*Dt
                    end do                            
                    DthetaL = Dtheta(1)*DistU(1) + Dtheta(2)*DistU(2) + Dtheta(3)*DistU(3)
                    !  twisting deform
                    do K = 1,3
                        DthetaT(K) = DthetaL*DistU(K)
                    end do
                    DthetaTL = sqrt(DthetaT(1)*DthetaT(1) + DthetaT(2)*DthetaT(2) + DthetaT(3)*DthetaT(3))
                    !  rolling deform
                    do K = 1,3
                        DthetaR(K) = Dtheta(K) - DthetaT(K)
                    end do
                    DthetaRL = sqrt(DthetaR(1)*DthetaR(1) + DthetaR(2)*DthetaR(2) + DthetaR(3)*DthetaR(3))

                    !  rolling moment of Particle J
                    do K = 1,3
                        rolling_moment(K) = Kr*DthetaR(K) + Mr(K)
                    end do
                    rolling_momentL = sqrt(rolling_moment(1)*rolling_moment(1) + rolling_moment(2)*rolling_moment(2) + rolling_moment(3)*rolling_moment(3))                

                    !if (rolling) then  !  Have rolled
                    !    if (DthetaRL .GT. 1.0e-14) then  !  Still slipping
                    !        do K = 1,3
                    !            rolling_moment(K) = 2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL  !  Particle J
                    !        end do
                    !    else  !  Approach sticking
                    !        do K = 1,3
                    !            rolling_moment(K) = 0.0D0  !  Particle J
                    !        end do
                    !        rolling = .false.
                    !    end if
                    !else
                        if (rolling_momentL .GT. 2.1D0*0.25D0*m_Beta*Rij*normal_forceL) then  !  Rolling
                            rolling = .true.
                            if (DthetaRL .GT. 1.0e-14) then
                                do K = 1,3
                                    rolling_moment(K) = 2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL
                                    rolling_history(K) = rolling_moment(K)
                                end do
                                !  Add frictional heat
                                Heat(2,particleI) = Heat(2,particleI) + 0.5D0*2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaRL
                                !Heat(2,particleJ) = Heat(2,particleJ) + 0.5D0*2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaRL
                            else
                                do K = 1,3
                                    rolling_moment(K) = 0.0D0
                                    rolling_history(K) = rolling_moment(K)
                                end do
                            end if
                        else
                            rolling = .false.  !  Sticking
                            do K = 1,3
                                rolling_history(K) = rolling_moment(K)
                                rolling_moment(K) = rolling_moment(K) - Cr*DthetaR(K)/Dt
                            end do
                        end if
                    !end if  

                    !  twisting moment of Particle J
                    do K = 1,3
                        twisting_moment(K) = Kt*DthetaT(K) + Mt(K)
                    end do
                    twisting_momentL = sqrt(twisting_moment(1)*twisting_moment(1) + twisting_moment(2)*twisting_moment(2) + twisting_moment(3)*twisting_moment(3))                

                    !if (twisting) then  !  Have twisted
                    !    if (DthetaTL .GT. 1.0e-14) then  !  Still slipping
                    !        do K = 1,3
                    !            twisting_moment(K) = 0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL  !  Particle J
                    !        end do
                    !    else  !  Approach sticking
                    !        do K = 1,3
                    !            twisting_moment(K) = 0.0D0  !  Particle J
                    !        end do
                    !        twisting = .false.
                    !    end if
                    !else
                        if (twisting_momentL .GT. 0.65D0*m_mu_s*m_Beta*Rij*normal_forceL) then  !  Rolling
                            twisting = .true.
                            if (DthetaTL .GT. 1.0e-14) then
                                do K = 1,3
                                    twisting_moment(K) = 0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL
                                    twisting_history(K) = twisting_moment(K)
                                end do
                                !  Add frictional heat
                                Heat(3,particleI) = Heat(3,particleI) + 0.5D0*0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaTL
                                !Heat(3,particleJ) = Heat(3,particleJ) + 0.5D0*0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaTL
                            else
                                do K = 1,3
                                    twisting_moment(K) = 0.0D0
                                    twisting_history(K) = twisting_moment(K)
                                end do
                            end if
                        else
                            twisting = .false.  !  Sticking
                            do K = 1,3
                                twisting_history(K) = twisting_moment(K)
                                twisting_moment(K) = twisting_moment(K) - Ct*DthetaT(K)/Dt
                            end do
                        end if
                    !end if  

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
                    !  Apply moment
                    do K = 1,3
                        !FM(K,particleJ) = + rolling_moment(K) + twisting_moment(K) + FM(K,particleJ) 
                        FM(K,particleI) = - rolling_moment(K) - twisting_moment(K) + FM(K,particleI) 
                    end do                                          
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

                    !if (particleI.EQ.PP .OR. particleJ.EQ.PP) then
                    !    write(123,'(2I5,2X,15F30.17)',advance='no') particleI,particleJ,(Vtot(K),K=1,3),DistS,DistL,ERR,(DistU(K),K=1,3),(Vnor(K),K=1,3),(Vtan(K),K=1,3)
                    !    write(123,'(18F30.17)') (rolling_moment(K),K=1,3),(DthetaR(K),K=1,3),(Mr(K),K=1,3),(twisting_moment(K),K=1,3),(DthetaT(K),K=1,3),(Mt(K),K=1,3)
                    !end if

                    !  cohesive force
                    do K = 1,3
                        cohesive_force(K) = - m_c*(m_Beta*Rij)**2*DistU(K)
                    end do
                    do K = 1,3
                        !F(K,particleJ) = F(K,particleJ) + cohesive_force(K)
                        F(K,particleI) = F(K,particleI) - cohesive_force(K)
                    end do
                    !  memory the contact in the Hertz linklist.
                    if (associated(Temp%prev)) then
                        !  Temp is in center of linklist!!!
                        if (Temp%No .EQ. particleJ) then
                            !  Have contacted.
                            do K = 1,3
                                Temp%Hertz(K) = tangential_history(K)
                                Temp%Mrot(K) = rolling_history(K)
                                Temp%Mtwist(K) = twisting_history(K)
                            end do
                            Temp%is_touching = touching
                            Temp%is_slipping = slipping
                            Temp%is_rolling = rolling
                            Temp%is_twisting = twisting
                            !Temp%recordTime = Time + Dt
                            Tail => Temp
                        else
                            !  First contacted.
                            allocate(TempH)
                            TempH = Nodelink(particleJ,tangential_history,rolling_history,twisting_history,&
                            & touching,slipping,rolling,twisting,Temp,Temp%next)
                            if (associated(Temp%next)) Temp%next%prev => TempH
                            Temp%next => TempH
                            Head(particleI)%No = LenNode + 1
                            Tail => TempH
                        end if
                    else
                        !  Temp is Head of linklist!!!
                        allocate(TempH)
                        TempH = Nodelink(particleJ,tangential_history,rolling_history,twisting_history,&
                        & touching,slipping,rolling,twisting,Temp,Temp%next)
                        if (associated(Temp%next)) Temp%next%prev => TempH
                        Temp%next => TempH
                        Head(particleI)%No = LenNode + 1
                        Tail => TempH
                    end if
                else
                    !  memory the separation in the Hertz linklist.
                    if (associated(Temp%prev) .AND. Temp%No.EQ.particleJ) then
                        !  Temp is center of linklist!!!
                        Temp%prev%next => Temp%next
                        if(associated(Temp%next)) Temp%next%prev => Temp%prev
                        Head(particleI)%No = LenNode - 1
                        Tail => Temp%prev
                        deallocate(Temp)
                        !  When else Temp is Head of linklist!!!
                    else
                        Tail => Temp
                    end if
                    Rij = R(particleI)*R(particleJ)/(R(particleI)+R(particleJ))
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                              
                    !  no contact but in region 1
                    if (Dn > -Rij*(m_r_cut-1.0D0)*0.5D0) then
                        DistR = 1.0D0/DistL
                        !  calculate the normal vector
                        do K=1,3
                            DistU(K) = Dist(K)*DistR
                        end do
                        !  cohesive force
                        do K = 1,3
                            cohesive_force(K) = - m_c*(m_Beta*Rij)**2*DistU(K)
                        end do
                        !  apply force
                        do K = 1,3
                            !F(K,particleJ) = F(K,particleJ) + cohesive_force(K)
                            F(K,particleI) = F(K,particleI) - cohesive_force(K)
                        end do
                    else if (Dn > -Rij*(m_r_cut-1.0D0)) then
                        DistR = 1.0D0/DistL
                        !  calculate the normal vector
                        do K=1,3
                            DistU(K) = Dist(K)*DistR
                        end do
                        !  cohesive force
                        do K = 1,3
                            cohesive_force(K) = - m_c*(m_Beta)**2*Rij &
                            &*2.0D0*(Dn/(m_r_cut-1.0D0) + Rij)*DistU(K)
                        end do
                        do K = 1,3
                            !F(K,particleJ) = F(K,particleJ) + cohesive_force(K)
                            F(K,particleI) = F(K,particleI) - cohesive_force(K)
                        end do 
                    end if
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                end if
            end if
            !if (particleI.EQ.PP .OR. particleJ.EQ.PP) then
            !    write(123,*)
            !end if
        end do     

        do JJ = 1,26
            Tail => Head(particleI)
            IDJ = DEM(Linklist(I))%NeighborID(JJ)
            if (IDJ .NE. 0) then
                IDOuterJ => IDInner(IDJ)
                do while(associated(IDOuterJ%next))
                    IDOuterJ => IDOuterJ%next
                    particleJ = IDOuterJ%No

                    !if (particleI.EQ.PP .OR. particleJ.EQ.PP) then
                    !    write(123,'(F15.5,2X,A5,2X,2I5,2X)',advance='no') Time,'outer',particleI,particleJ
                    !end if

                    !if (particleI.EQ.PP .OR. particleJ.EQ.PP) then
                    !    write(123,'(A6,2X,I5,2X)',advance='no') 'check:',check
                    !end if

                    if (.true.) then   !  .OR. J.EQ.N
                        !  Initialize state params
                        do K = 1,3
                            Dist(K) = X(K,particleJ) -  X(K,particleI)
                            H(K) = 0.0D0
                            Mr(K) = 0.0D0
                            Mt(K) = 0.0D0
                        end do
                        touching = .false.
                        slipping = .false.
                        rolling = .false.
                        twisting = .false.
                        !  Distance vector
                        DistS = Dist(1)*Dist(1) + Dist(2)*Dist(2) + Dist(3)*Dist(3)
                        DistL = sqrt(DistS)
                        Dn = R(particleI) + R(particleJ) - DistL
                        !  lookup this contact(or not) in the Hertz list I
                        LenNode = Head(particleI)%No
                        Temp => Tail
                        do while(associated(Temp%next))
                            Temp => Temp%next
                            if (Temp%No .EQ. particleJ) then
                                do K = 1,3
                                    H(K) = Temp%Hertz(K)
                                    Mr(K) = Temp%Mrot(K)
                                    Mt(K) = Temp%Mtwist(K)
                                end do
                                touching = Temp%is_touching
                                slipping = Temp%is_slipping
                                rolling = Temp%is_rolling
                                twisting = Temp%is_twisting
                                exit
                            else if (Temp%No .GT. particleJ) then
                                Temp => Temp%prev
                                exit
                            end if
                        end do        
                        !if (particleI.EQ.PP .OR. particleJ.EQ.PP) then
                        !    write(123,'(F15.5,2X)',advance='no') Dn
                        !end if
                        !  When collision calculate the repulsive restoring spring force which is generated along the normal and tangential according to Hooke's law
                        if (Dn .GT. 0.0D0) then
                            DistR = 1.0D0/DistL
                            !  calculate the normal vector
                            do K=1,3
                                DistU(K) = Dist(K)*DistR
                            end do
                            Ap = (R(particleI)*R(particleI)-R(particleJ)*R(particleJ)+DistS)/2.0D0*DistR
                            An = DistL-Ap
#ifdef HertzMindlinVisco                    
                            !  calculate material constant
                            Rij = R(particleI)*R(particleJ)/(R(particleI)+R(particleJ))
                            Mij = Body(particleI)*Body(particleJ)/(Body(particleI)+Body(particleJ))
                            Kn = 2.0D0*m_E*sqrt(Rij*Dn)/(3.0D0*(1.0D0-m_nu*m_nu))
                            Cn = -Kn*m_A
                            Ks = 2.0D0*m_E/(1.0D0+m_nu)/(2.0D0-m_nu)*sqrt(Rij)*sqrt(Dn)
                            !  select tangential damping mode
                            if (m_COR > 1.0D0) then
                                Cs = -2.0D0*m_E/(1.0D0+m_nu)/(2.0D0-m_nu)*sqrt(Dn)*m_A
                            elseif (m_COR >= 0.0D0) then
                                lnCOR=log(m_COR)
                                Cs = 2.0D0*sqrt(5.0D0/6.0D0)*lnCOR/sqrt(lnCOR**2+3.1415926D0**2) &
                                & *sqrt(2.0D0*Mij*m_E/(1.0D0+m_nu)/(2.0D0-m_nu))*(Rij**0.25)*(Dn**0.25)
                            else
                                Cs = 0.0D0
                            end if
                            Kr = 0.25D0*Kn*(m_Beta*Rij)**2
                            Cr = 0.25D0*Cn*(m_Beta*Rij)**2
                            Kt = 0.5D0*Ks*(m_Beta*Rij)**2
                            Ct = 0.5D0*Cs*(m_Beta*Rij)**2
#elif HertzMindlinResti
                            !  calculate material constant
                            Rij = R(particleI)*R(particleJ)/(R(particleI)+R(particleJ))
                            Mij = Body(particleI)*Body(particleJ)/(Body(particleI)+Body(particleJ))
                            Kn = 2.0D0*m_E*sqrt(Rij*Dn)/(3.0D0*(1.0D0-m_nu*m_nu))
                            lnCOR = log(m_COR)
                            Cn = 2.0D0*sqrt(5.0D0/6.0D0)*lnCOR/sqrt(lnCOR**2+3.1415926D0**2) &
                            & *sqrt(Mij*m_E/(1.0D0-m_nu*m_nu))*(Rij**0.25)*(Dn**0.25)
                            Ks = 2.0D0*m_E/(1.0D0+m_nu)/(2.0D0-m_nu)*sqrt(Rij)*sqrt(Dn)
                            Cs = 2.0D0*sqrt(5.0D0/6.0D0)*lnCOR/sqrt(lnCOR**2+3.1415926D0**2) &
                            & *sqrt(2.0D0*Mij*m_E/(1.0D0+m_nu)/(2.0D0-m_nu))*(Rij**0.25)*(Dn**0.25)
                            Kr = 0.25D0*Kn*(m_Beta*Rij)**2
                            Cr = 0.25D0*Cn*(m_Beta*Rij)**2
                            Kt = 0.5D0*Ks*(m_Beta*Rij)**2
                            Ct = 0.5D0*Cs*(m_Beta*Rij)**2
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
                                normal_force(K) = Kn*Dn*DistU(K) + Cn*Vnor(K)
                            end do
                            normal_forceL = sqrt(normal_force(1)*normal_force(1) + normal_force(2)*normal_force(2) + normal_force(3)*normal_force(3))

                            !  Add energy
                            Energy(particleI) = Energy(particleI) + 0.2D0*Kn*(Dn**2)
                            !Energy(particleJ) = Energy(particleJ) + 0.2D0*Kn*(Dn**2)

                            !  tangential deform
                            do K = 1,3
                                Ds(K) = Vtan(K)*Dt
                            end do
                            DsL = sqrt(Ds(1)*Ds(1) + Ds(2)*Ds(2) + Ds(3)*Ds(3))

                            !  tangential force of Particle J
                            do K = 1,3
                                tangential_force(K) = - Ks*Ds(K) + Cs*Vtan(K) + H(K)
                            end do
                            tangential_forceL = sqrt(tangential_force(1)*tangential_force(1) + tangential_force(2)*tangential_force(2) + tangential_force(3)*tangential_force(3))

                            !if (particleI.EQ.PP .OR. particleJ.EQ.PP) then
                            !    write(123,'(14F15.5,2X)',advance='no') normal_forceL,tangential_forceL,(H(K),K=1,3),(Vtan(K),K=1,3),(Ds(K),K=1,3),(tangential_force(K),K=1,3)
                            !end if

                            !if (slipping) then  !  Have slipped
                            !    if (DsL .GT. 1.0e-14) then  !  Still slipping
                            !        do K = 1,3
                            !            tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL  !  Particle J
                            !        end do
                            !    else  !  Approach sticking
                            !        do K = 1,3
                            !            tangential_force(K) = 0.0D0  !  Particle J
                            !        end do
                            !        slipping = .false.
                            !    end if
                            !else
                                if (tangential_forceL .GT. normal_forceL*m_mu_s) then  !  Slipping
                                    slipping = .true.
                                    if (DsL .GT. 1.0e-14) then
                                        do K = 1,3
                                            tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL
                                            tangential_history(K) = tangential_force(K)
                                        end do
                                        !  Add frictional heat
                                        Heat(1,particleI) = Heat(1,particleI) + 0.5D0*m_mu_d*normal_forceL*DsL
                                        !Heat(1,particleJ) = Heat(1,particleJ) + 0.5D0*m_mu_d*normal_forceL*DsL
                                    else
                                        do K = 1,3
                                            tangential_force(K) = 0.0D0
                                            tangential_history(K) = tangential_force(K)
                                        end do
                                    end if
                                else
                                    slipping = .false.  !  Sticking
                                    do K = 1,3
                                        tangential_history(K) = tangential_force(K) - Cs*Vtan(K)
                                    end do
                                end if
                            !end if

                            touching = .true.
                            !  Apply force
                            do K = 1,3
                                !F(K,particleJ) = + normal_force(K) + tangential_force(K) + F(K,particleJ)
                                F(K,particleI) = - normal_force(K) - tangential_force(K) + F(K,particleI)
                            end do
                            !  Apply moment
                            !FM(1,particleJ) = - An*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,particleJ) 
                            !FM(2,particleJ) = - An*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,particleJ) 
                            !FM(3,particleJ) = - An*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,particleJ)                     
                            FM(1,particleI) = - Ap*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,particleI) 
                            FM(2,particleI) = - Ap*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,particleI) 
                            FM(3,particleI) = - Ap*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,particleI)

                            !  rolling
                            do K = 1,3
                                Dtheta(K) = (W(K,particleI)-W(K,particleJ))*Dt
                            end do                            
                            DthetaL = Dtheta(1)*DistU(1) + Dtheta(2)*DistU(2) + Dtheta(3)*DistU(3)
                            !  twisting deform
                            do K = 1,3
                                DthetaT(K) = DthetaL*DistU(K)
                            end do
                            DthetaTL = sqrt(DthetaT(1)*DthetaT(1) + DthetaT(2)*DthetaT(2) + DthetaT(3)*DthetaT(3))
                            !  rolling deform
                            do K = 1,3
                                DthetaR(K) = Dtheta(K) - DthetaT(K)
                            end do
                            DthetaRL = sqrt(DthetaR(1)*DthetaR(1) + DthetaR(2)*DthetaR(2) + DthetaR(3)*DthetaR(3))

                            !  rolling moment of Particle J
                            do K = 1,3
                                rolling_moment(K) = Kr*DthetaR(K) + Mr(K)
                            end do
                            rolling_momentL = sqrt(rolling_moment(1)*rolling_moment(1) + rolling_moment(2)*rolling_moment(2) + rolling_moment(3)*rolling_moment(3))                

                            !if (rolling) then  !  Have rolled
                            !    if (DthetaRL .GT. 1.0e-14) then  !  Still slipping
                            !        do K = 1,3
                            !            rolling_moment(K) = 2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL  !  Particle J
                            !        end do
                            !    else  !  Approach sticking
                            !        do K = 1,3
                            !            rolling_moment(K) = 0.0D0  !  Particle J
                            !        end do
                            !        rolling = .false.
                            !    end if
                            !else
                                if (rolling_momentL .GT. 2.1D0*0.25D0*m_Beta*Rij*normal_forceL) then  !  Rolling
                                    rolling = .true.
                                    if (DthetaRL .GT. 1.0e-14) then
                                        do K = 1,3
                                            rolling_moment(K) = 2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL
                                            rolling_history(K) = rolling_moment(K)
                                        end do
                                        !  Add frictional heat
                                        Heat(2,particleI) = Heat(2,particleI) + 0.5D0*2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaRL
                                        !Heat(2,particleJ) = Heat(2,particleJ) + 0.5D0*2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaRL
                                    else
                                        do K = 1,3
                                            rolling_moment(K) = 0.0D0
                                            rolling_history(K) = rolling_moment(K)
                                        end do
                                    end if
                                else
                                    rolling = .false.  !  Sticking
                                    do K = 1,3
                                        rolling_history(K) = rolling_moment(K)
                                        rolling_moment(K) = rolling_moment(K) - Cr*DthetaR(K)/Dt
                                    end do
                                end if
                            !end if     

                            !  twisting moment of Particle J
                            do K = 1,3
                                twisting_moment(K) = Kt*DthetaT(K) + Mt(K)
                            end do
                            twisting_momentL = sqrt(twisting_moment(1)*twisting_moment(1) + twisting_moment(2)*twisting_moment(2) + twisting_moment(3)*twisting_moment(3))                

                            !if (twisting) then  !  Have twisted
                            !    if (DthetaTL .GT. 1.0e-14) then  !  Still slipping
                            !        do K = 1,3
                            !            twisting_moment(K) = 0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL  !  Particle J
                            !        end do
                            !    else  !  Approach sticking
                            !        do K = 1,3
                            !            twisting_moment(K) = 0.0D0  !  Particle J
                            !        end do
                            !        twisting = .false.
                            !    end if
                            !else
                                if (twisting_momentL .GT. 0.65D0*m_mu_s*m_Beta*Rij*normal_forceL) then  !  Rolling
                                    twisting = .true.
                                    if (DthetaTL .GT. 1.0e-14) then
                                        do K = 1,3
                                            twisting_moment(K) = 0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL
                                            twisting_history(K) = twisting_moment(K)
                                        end do
                                        !  Add frictional heat
                                        Heat(3,particleI) = Heat(3,particleI) + 0.5D0*0.65D0*m_mu_s*m_Beta*Rij*normal_forceL*DthetaTL
                                        !Heat(3,particleJ) = Heat(3,particleJ) + 0.5D0*0.65D0*m_mu_s*m_Beta*Rij*normal_forceL*DthetaTL
                                    else
                                        do K = 1,3
                                            twisting_moment(K) = 0.0D0
                                            twisting_history(K) = twisting_moment(K)
                                        end do
                                    end if
                                else
                                    twisting = .false.  !  Sticking
                                    do K = 1,3
                                        twisting_history(K) = twisting_moment(K)
                                        twisting_moment(K) = twisting_moment(K) - Ct*DthetaT(K)/Dt
                                    end do
                                end if
                            !end if  

                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
                            !  Apply moment
                            do K = 1,3
                                !FM(K,particleJ) = + rolling_moment(K) + twisting_moment(K) + FM(K,particleJ) 
                                FM(K,particleI) = - rolling_moment(K) - twisting_moment(K) + FM(K,particleI) 
                            end do                                          
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

                            !if (particleI.EQ.PP .OR. particleJ.EQ.PP) then
                            !    write(123,'(2I5,2X,15F30.17)',advance='no') particleI,particleJ,(Vtot(K),K=1,3),DistS,DistL,ERR,(DistU(K),K=1,3),(Vnor(K),K=1,3),(Vtan(K),K=1,3)
                            !    write(123,'(18F30.17)') (rolling_moment(K),K=1,3),(DthetaR(K),K=1,3),(Mr(K),K=1,3),(twisting_moment(K),K=1,3),(DthetaT(K),K=1,3),(Mt(K),K=1,3)
                            !end if

                            !  cohesive force
                            do K = 1,3
                                cohesive_force(K) = - m_c*(m_Beta*Rij)**2*DistU(K)
                            end do
                            do K = 1,3
                                !F(K,particleJ) = F(K,particleJ) + cohesive_force(K)
                                F(K,particleI) = F(K,particleI) - cohesive_force(K)
                            end do
                            !  memory the contact in the Hertz linklist.
                            if (associated(Temp%prev)) then
                                !  Temp is in center of linklist!!!
                                if (Temp%No .EQ. particleJ) then
                                    !  Have contacted.
                                    do K = 1,3
                                        Temp%Hertz(K) = tangential_history(K)
                                        Temp%Mrot(K) = rolling_history(K)
                                        Temp%Mtwist(K) = twisting_history(K)
                                    end do
                                    Temp%is_touching = touching
                                    Temp%is_slipping = slipping
                                    Temp%is_rolling = rolling
                                    Temp%is_twisting = twisting
                                    !Temp%recordTime = Time + Dt
                                    Tail => Temp
                                else
                                    !  First contacted.
                                    allocate(TempH)
                                    TempH = Nodelink(particleJ,tangential_history,rolling_history,twisting_history,&
                                    & touching,slipping,rolling,twisting,Temp,Temp%next)
                                    if (associated(Temp%next)) Temp%next%prev => TempH
                                    Temp%next => TempH
                                    Head(particleI)%No = LenNode + 1
                                    Tail => TempH
                                end if
                            else
                                !  Temp is Head of linklist!!!
                                allocate(TempH)
                                TempH = Nodelink(particleJ,tangential_history,rolling_history,twisting_history,&
                                & touching,slipping,rolling,twisting,Temp,Temp%next)
                                if (associated(Temp%next)) Temp%next%prev => TempH
                                Temp%next => TempH
                                Head(particleI)%No = LenNode + 1
                                Tail => TempH
                            end if
                        else
                            !  memory the separation in the Hertz linklist.
                            if (associated(Temp%prev) .AND. Temp%No.EQ.particleJ) then
                                !  Temp is center of linklist!!!
                                Temp%prev%next => Temp%next
                                if(associated(Temp%next)) Temp%next%prev => Temp%prev
                                Head(particleI)%No = LenNode - 1
                                Tail => Temp%prev
                                deallocate(Temp)
                                !  When else Temp is Head of linklist!!!
                            else
                                Tail => Temp
                            end if
                            Rij = R(particleI)*R(particleJ)/(R(particleI)+R(particleJ))
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            if (Dn > -Rij*(m_r_cut-1.0D0)*0.5D0) then
                                DistR = 1.0D0/DistL
                                !  calculate the normal vector
                                do K=1,3
                                    DistU(K) = Dist(K)*DistR
                                end do
                                !  cohesive force
                                do K = 1,3
                                    cohesive_force(K) = - m_c*(m_Beta*Rij)**2*DistU(K)
                                end do
                                !  apply force
                                do K = 1,3
                                    !F(K,particleJ) = F(K,particleJ) + cohesive_force(K)
                                    F(K,particleI) = F(K,particleI) - cohesive_force(K)
                                end do                                     
                            else if (Dn > -Rij*(m_r_cut-1.0D0)) then
                                DistR = 1.0D0/DistL
                                !  calculate the normal vector
                                do K=1,3
                                    DistU(K) = Dist(K)*DistR
                                end do
                                !  cohesive force
                                do K = 1,3
                                    cohesive_force(K) = - m_c*m_Beta**2*Rij &
                                    &*2.0D0*(Dn/(m_r_cut-1.0D0) + Rij)*DistU(K)
                                end do
                                do K = 1,3
                                    !F(K,particleJ) = F(K,particleJ) + cohesive_force(K)
                                    F(K,particleI) = F(K,particleI) - cohesive_force(K)
                                end do      
                            end if
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        end if
                    end if
                    !if (particleI.EQ.PP .OR. particleJ.EQ.PP) then
                    !    write(123,*)
                    !end if
                end do  
            end if
        end do     
    end do
    !$OMP END PARALLEL DO
    !oend = omp_get_wtime()
    !write(*,*) "forceParticleLattice",(oend-ostart)
    
    !write(FileNameForce,'(F10.5)') Time
    !FileNameForce = trim(FileNameForce)//'Force.txt'
    !open(10000,FILE=FileNameForce)
    !do I = 1,N
    !    write(10000,'(3F30.15)') (F(K,I),K=1,3)
    !end do
    !close(10000)
    !close(123)
    !pause
    !
    !stop

    return
    end
