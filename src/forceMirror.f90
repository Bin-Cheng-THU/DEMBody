    !********************************************************************
    !     DEMBody 4.4
    !     ***********
    !
    !     Force for mirrored particles.
    !     --------------------------
    !      
    !     @Using rolling friction similar to LIGGGHTS
    !     @Using Hertz-Mindlin contact model similar to Wada
    !     @Using Cohesion model similar to Scheeres
    !     @Using damping model applicable for ice ball
    !
    !********************************************************************
    subroutine forceMirror()

    use global
    implicit none

    real(kind=8)  Dist(3),DistS,DistL,DistR,DistU(3)
    real(kind=8)  Vrel(3),Vrot(3),Vtot(3),ERR,Vnor(3),Vtan(3)
    real(kind=8)  normal_force(3),normal_forceL
    real(kind=8)  tangential_force(3),tangential_forceL
    real(kind=8)  rolling_moment(3),rolling_momentL
    real(kind=8)  twisting_moment(3),twisting_momentL
    real(kind=8)  cohesive_force(3)
    real(kind=8)  Ap,An
    real(kind=8)  Rij,Mij,Iij
    real(kind=8)  Kn,Cn,Ks,Cs,Kr,Cr,Kt,Ct,lnCOR
    real(kind=8)  Dn,Ds(3),DsL,Dtheta(3),DthetaL,DthetaR(3),DthetaRL,DthetaT(3),DthetaTL
    real(kind=8)  H(3),Mr(3),Mt(3)
    real(kind=8)  RV(3)
    logical :: slipping,rolling,twisting     !  State of friction of T
    logical :: touching                      !  State of touch of T
    integer :: I,J,K,L,LenNode               !  Iterator
    type(Nodelink),pointer :: Temp           !  Temporary pointer
    type(Nodelink),pointer :: TempH          !  Contact pointer
    type(Nodelink),pointer :: Tail           !  Tail pointer
    real(kind=8) ShearPBC                    !  shearPBC                    

    !character(30) :: FileName
    !write(FileName,'(F10.5)') Time
    !FileName = trim(FileName)//'HeadMirror.txt'
    !open(123,FILE=FileName)
    
    !  Loop over all bodies through the NodeTree.
    !$OMP PARALLEL DO &
    !$OMP& PRIVATE(shearPBC,Temp,TempH,Tail,LenNode,&
    !$OMP& I,J,K,Dist,DistS,DistL,DistR,DistU,Vrel,Vrot,Vtot,ERR,Vnor,Vtan,&
    !$OMP& normal_force,normal_forceL,tangential_force,tangential_forceL,&
    !$OMP& rolling_moment,rolling_momentL,twisting_moment,twisting_momentL,cohesive_force,Ap,An,Rij,Mij,Iij,&
    !$OMP& Kn,Cn,Ks,Cs,Kr,Cr,Kt,Ct,lnCOR,Dn,Ds,DsL,Dtheta,DthetaL,DthetaR,DthetaRL,DthetaT,DthetaTL,H,Mr,Mt,RV,&
    !$OMP& slipping,rolling,twisting,touching) SCHEDULE(GUIDED)
    do I = 1,N
        
        !################         Wall X1          ################### 
        if (Tag1(I).EQ.1) then
            Tail => Head(I)
            do J = 1,N
                if (Tag2(J).EQ.1) then
                    
                    !if (I.EQ.PP .OR. J.EQ.PP) then
                    !    write(123,'(F15.5,2X,2I5,2X)',advance='no') Time,I,J
                    !end if
            
                    !  Initialize state params
                    do K = 1,3
                        Dist(K) = X(K,J) -  X(K,I)
                        H(K) = 0.0D0
                        Mr(K) = 0.0D0
                        Mt(K) = 0.0D0
                    end do
                    Dist(1) = Dist(1) - LenBoxX
                    touching = .false.
                    slipping = .false.
                    rolling = .false.
                    twisting = .false.
                    !  Distance vector
                    DistS = Dist(1)*Dist(1) + Dist(2)*Dist(2) + Dist(3)*Dist(3)
                    DistL = sqrt(DistS)
                    Dn = R(I) + R(J) - DistL
                    !  lookup this contact(or not) in the Hertz list
                    LenNode = Head(I)%No
                    Temp => Tail
                    do while(associated(Temp%next))
                        Temp => Temp%next
                        if (Temp%No .EQ. J) then
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
                        else if (Temp%No .GT. J) then
                            Temp => Temp%prev
                            exit
                        end if
                    end do
                    !if (I.EQ.PP .OR. J.EQ.PP) then
                    !    write(123,'(F15.5,2X)',advance='no') Dn
                    !end if
                    !  When collision calculate the repulsive restoring spring force which is generated along the normal and tangential according to Hooke's law
                    if (Dn .GT. 0.0D0) then
                        DistR = 1.0D0/DistL
                        !  calculate the normal vector
                        do K=1,3
                            DistU(K) = Dist(K)*DistR
                        end do
                        Ap = (R(I)*R(I)-R(J)*R(J)+DistS)/2.0D0*DistR
                        An = DistL-Ap
#ifdef HertzMindlinVisco                    
                        !  calculate material constant
                        Rij = R(I)*R(J)/(R(I)+R(J))
                        Mij = Body(I)*Body(J)/(Body(I)+Body(J))
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
                        Rij = R(I)*R(J)/(R(I)+R(J))
                        Mij = Body(I)*Body(J)/(Body(I)+Body(J))
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
                            Vrel(K) = Xdot(K,J) - Xdot(K,I)
                        end do
                        !  negative rotate relative velocity
                        Vrot(1) = (DistU(2)*W(3,J) - DistU(3)*W(2,J))*An + (DistU(2)*W(3,I) - DistU(3)*W(2,I))*Ap
                        Vrot(2) = (DistU(3)*W(1,J) - DistU(1)*W(3,J))*An + (DistU(3)*W(1,I) - DistU(1)*W(3,I))*Ap
                        Vrot(3) = (DistU(1)*W(2,J) - DistU(2)*W(1,J))*An + (DistU(1)*W(2,I) - DistU(2)*W(1,I))*Ap
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
                        Energy(I) = Energy(I) + 0.4D0*Kn*(Dn**2)

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

                        !if (I.EQ.PP .OR. J.EQ.PP) then
                        !    write(123,'(11F15.5,2X)',advance='no') normal_forceL,tangential_forceL,(H(K),K=1,3),(Vtan(K),K=1,3),(Ds(K),K=1,3)
                        !end if
                
                        if (slipping) then  !  Have slipped
                            if (DsL .GT. 1.0e-8) then  !  Still slipping
                                do K = 1,3
                                    tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    tangential_force(K) = 0.0D0  !  Particle J
                                end do
                                slipping = .false.
                            end if
                        else
                            if (tangential_forceL .GT. normal_forceL*m_mu_s) then  !  Slipping
                                slipping = .true.
                                if (DsL .GT. 1.0e-14) then
                                    do K = 1,3
                                        tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL
                                    end do
                                else
                                    do K = 1,3
                                        tangential_force(K) = 0.0D0
                                    end do
                                end if
                            else
                                slipping = .false.  !  Sticking
                            end if
                        end if

                        touching = .true.
                        !  Apply force
                        do K = 1,3
                            F(K,I) = - normal_force(K) - tangential_force(K) + F(K,I)
                        end do
                        !  Apply moment                     
                        FM(1,I) = - Ap*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,I) 
                        FM(2,I) = - Ap*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,I) 
                        FM(3,I) = - Ap*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,I)

                        !  rolling
                        do K = 1,3
                            Dtheta(K) = (W(K,I)-W(K,J))*Dt
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

                        if (rolling) then  !  Have rolled
                            if (DthetaRL .GT. 1.0e-8) then  ! Still slipping
                                do K = 1,3
                                    rolling_moment(K) = 2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    rolling_moment(K) = 0.0D0  !  Particle J
                                end do
                                rolling = .false.
                            end if
                        else
                            if (rolling_momentL .GT. 2.1D0*0.25D0*m_Beta*Rij*normal_forceL) then  !  Rolling
                                rolling = .true.
                                if (DthetaRL .GT. 1.0e-14) then
                                    do K = 1,3
                                        rolling_moment(K) = 2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL
                                    end do
                                else
                                    do K = 1,3
                                        rolling_moment(K) = 0.0D0
                                    end do
                                end if
                            else
                                rolling = .false.  !  Sticking
                                do K = 1,3
                                    rolling_moment(K) = rolling_moment(K) - Cr*DthetaR(K)/Dt
                                end do
                            end if
                        end if                                  

                        !  twisting moment of Particle J
                        do K = 1,3
                            twisting_moment(K) = Kt*DthetaT(K) + Mt(K)
                        end do
                        twisting_momentL = sqrt(twisting_moment(1)*twisting_moment(1) + twisting_moment(2)*twisting_moment(2) + twisting_moment(3)*twisting_moment(3))                

                        if (twisting) then  !  Have twisted
                            if (DthetaTL .GT. 1.0e-8) then  !  Still slipping
                                do K = 1,3
                                    twisting_moment(K) = 0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    twisting_moment(K) = 0.0D0  !  Particle J
                                end do
                                twisting = .false.
                            end if
                        else
                            if (twisting_momentL .GT. 0.65D0*m_mu_s*m_Beta*Rij*normal_forceL) then  !  Rolling
                                twisting = .true.
                                if (DthetaTL .GT. 1.0e-14) then
                                    do K = 1,3
                                        twisting_moment(K) = 0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL
                                    end do
                                else
                                    do K = 1,3
                                        twisting_moment(K) = 0.0D0
                                    end do
                                end if
                            else
                                twisting = .false.  !  Sticking
                                do K = 1,3
                                    twisting_moment(K) = twisting_moment(K) - Ct*DthetaT(K)/Dt
                                end do
                            end if
                        end if  

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
                        !  Apply moment
                        do K = 1,3
                            FM(K,I) = - rolling_moment(K) - twisting_moment(K) + FM(K,I) 
                        end do
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        
                        !  cohesive force
                        do K = 1,3
                            cohesive_force(K) = - m_c*(m_Beta*Rij)**2*DistU(K)
                        end do
                        do K = 1,3
                            F(K,I) = F(K,I) - cohesive_force(K)
                        end do
                        !  memory the contact in the Hertz linklist.
                        if (associated(Temp%prev)) then
                            !  Temp is in center of linklist!!!
                            if (Temp%No .EQ. J) then
                                !  Have contacted.
                                do K = 1,3
                                    Temp%Hertz(K) = tangential_force(K)
                                    Temp%Mrot(K) = rolling_moment(K)
                                    Temp%Mtwist(K) = twisting_moment(K)
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
                                TempH = Nodelink(J,tangential_force,rolling_moment,twisting_moment,&
                                & touching,slipping,rolling,twisting,Temp,Temp%next)
                                if (associated(Temp%next)) Temp%next%prev => TempH
                                Temp%next => TempH
                                Head(I)%No = LenNode + 1
                                Tail => TempH
                            end if
                        else
                            !  Temp is Head of linklist!!!
                            allocate(TempH)
                            TempH = Nodelink(J,tangential_force,rolling_moment,twisting_moment,&
                            & touching,slipping,rolling,twisting,Temp,Temp%next)
                            if (associated(Temp%next)) Temp%next%prev => TempH
                            Temp%next => TempH
                            Head(I)%No = LenNode + 1
                            Tail => TempH
                        end if
                    else
                        !  memory the separation in the Hertz linklist.
                        if (associated(Temp%prev) .AND. Temp%No.EQ.J) then
                            !  Temp is center of linklist!!!
                            Temp%prev%next => Temp%next
                            if(associated(Temp%next)) Temp%next%prev => Temp%prev
                            Head(I)%No = LenNode - 1
                            Tail => Temp%prev
                            deallocate(Temp)
                            !  When else Temp is Head of linklist!!!
                        else
                            Tail => Temp                            
                        end if
                        Rij = R(I)*R(J)/(R(I)+R(J))
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
                                F(K,I) = F(K,I) - cohesive_force(K)
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
                            !  apply force
                            do K = 1,3
                                F(K,I) = F(K,I) - cohesive_force(K)
                            end do                                   
                        end if
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    end if
                    !if (I.EQ.PP .OR. J.EQ.PP) then
                    !    write(123,*)
                    !end if
                end if
            end do
        end if

        !################         Wall X2          ################### 
        if (Tag2(I).EQ.1) then
            Tail => Head(I)
            do J = 1,N
                if (Tag1(J).EQ.1) then   !  .OR. I.EQ.N
                    
                    !if (I.EQ.PP .OR. J.EQ.PP) then
                    !    write(123,'(F15.5,2X,2I5,2X)',advance='no') Time,I,J
                    !end if
                    
                    !  Initialize state params
                    do K = 1,3
                        Dist(K) = X(K,J) -  X(K,I)
                        H(K) = 0.0D0
                        Mr(K) = 0.0D0
                        Mt(K) = 0.0D0
                    end do
                    Dist(1) = Dist(1) + LenBoxX
                    touching = .false.
                    slipping = .false.
                    rolling = .false.
                    twisting = .false.
                    !  Distance vector
                    DistS = Dist(1)*Dist(1) + Dist(2)*Dist(2) + Dist(3)*Dist(3)
                    DistL = sqrt(DistS)
                    Dn = R(I) + R(J) - DistL
                    !  lookup this contact(or not) in the Hertz list
                    LenNode = Head(I)%No
                    Temp => Tail
                    do while(associated(Temp%next))
                        Temp => Temp%next
                        if (Temp%No .EQ. J) then
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
                        else if (Temp%No .GT. J) then
                            Temp => Temp%prev
                            exit
                        end if
                    end do
                    !if (I.EQ.PP .OR. J.EQ.PP) then
                    !    write(123,'(F15.5,2X)',advance='no') Dn
                    !end if
                    !  When collision calculate the repulsive restoring spring force which is generated along the normal and tangential according to Hooke's law
                    if (Dn .GT. 0.0D0) then
                        DistR = 1.0D0/DistL
                        !  calculate the normal vector
                        do K=1,3
                            DistU(K) = Dist(K)*DistR
                        end do
                        Ap = (R(I)*R(I)-R(J)*R(J)+DistS)/2.0D0*DistR
                        An = DistL-Ap
#ifdef HertzMindlinVisco                    
                        !  calculate material constant
                        Rij = R(I)*R(J)/(R(I)+R(J))
                        Mij = Body(I)*Body(J)/(Body(I)+Body(J))
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
                        Rij = R(I)*R(J)/(R(I)+R(J))
                        Mij = Body(I)*Body(J)/(Body(I)+Body(J))
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
                            Vrel(K) = Xdot(K,J) - Xdot(K,I)
                        end do
                        !  negative rotate relative velocity
                        Vrot(1) = (DistU(2)*W(3,J) - DistU(3)*W(2,J))*An + (DistU(2)*W(3,I) - DistU(3)*W(2,I))*Ap
                        Vrot(2) = (DistU(3)*W(1,J) - DistU(1)*W(3,J))*An + (DistU(3)*W(1,I) - DistU(1)*W(3,I))*Ap
                        Vrot(3) = (DistU(1)*W(2,J) - DistU(2)*W(1,J))*An + (DistU(1)*W(2,I) - DistU(2)*W(1,I))*Ap
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
                        Energy(I) = Energy(I) + 0.4D0*Kn*(Dn**2)

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

                        !if (I.EQ.PP .OR. J.EQ.PP) then
                        !    write(123,'(11F15.5,2X)',advance='no') normal_forceL,tangential_forceL,(H(K),K=1,3),(Vtan(K),K=1,3),(Ds(K),K=1,3)
                        !end if
                        
                        if (slipping) then  !  Have slipped
                            if (DsL .GT. 1.0e-8) then  !  Still slipping
                                do K = 1,3
                                    tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    tangential_force(K) = 0.0D0  !  Particle J
                                end do
                                slipping = .false.
                            end if
                        else
                            if (tangential_forceL .GT. normal_forceL*m_mu_s) then  !  Slipping
                                slipping = .true.
                                if (DsL .GT. 1.0e-14) then
                                    do K = 1,3
                                        tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL
                                    end do
                                else
                                    do K = 1,3
                                        tangential_force(K) = 0.0D0
                                    end do
                                end if
                            else
                                slipping = .false.  !  Sticking
                            end if
                        end if

                        touching = .true.
                        !  Apply force
                        do K = 1,3
                            F(K,I) = - normal_force(K) - tangential_force(K) + F(K,I)
                        end do
                        !  Apply moment                     
                        FM(1,I) = - Ap*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,I) 
                        FM(2,I) = - Ap*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,I) 
                        FM(3,I) = - Ap*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,I)

                        !  rolling
                        do K = 1,3
                            Dtheta(K) = (W(K,I)-W(K,J))*Dt
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

                        if (rolling) then  !  Have rolled
                            if (DthetaRL .GT. 1.0e-8) then  !  Still slipping
                                do K = 1,3
                                    rolling_moment(K) = 2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    rolling_moment(K) = 0.0D0  !  Particle J
                                end do
                                rolling = .false.
                            end if
                        else
                            if (rolling_momentL .GT. 2.1D0*0.25D0*m_Beta*Rij*normal_forceL) then  !  Rolling
                                rolling = .true.
                                if (DthetaRL .GT. 1.0e-14) then
                                    do K = 1,3
                                        rolling_moment(K) = 2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL
                                    end do
                                else
                                    do K = 1,3
                                        rolling_moment(K) = 0.0D0
                                    end do
                                end if
                            else
                                rolling = .false.  !  Sticking
                                do K = 1,3
                                    rolling_moment(K) = rolling_moment(K) - Cr*DthetaR(K)/Dt
                                end do
                            end if
                        end if                                  

                        !  twisting moment of Particle J
                        do K = 1,3
                            twisting_moment(K) = Kt*DthetaT(K) + Mt(K)
                        end do
                        twisting_momentL = sqrt(twisting_moment(1)*twisting_moment(1) + twisting_moment(2)*twisting_moment(2) + twisting_moment(3)*twisting_moment(3))                

                        if (twisting) then  !  Have twisted
                            if (DthetaTL .GT. 1.0e-8) then  !  Still slipping
                                do K = 1,3
                                    twisting_moment(K) = 0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    twisting_moment(K) = 0.0D0  !  Particle J
                                end do
                                twisting = .false.
                            end if
                        else
                            if (twisting_momentL .GT. 0.65D0*m_mu_s*m_Beta*Rij*normal_forceL) then  !  Rolling
                                twisting = .true.
                                if (DthetaTL .GT. 1.0e-14) then
                                    do K = 1,3
                                        twisting_moment(K) = 0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL
                                    end do
                                else
                                    do K = 1,3
                                        twisting_moment(K) = 0.0D0
                                    end do
                                end if
                            else
                                twisting = .false.  !  Sticking
                                do K = 1,3
                                    twisting_moment(K) = twisting_moment(K) - Ct*DthetaT(K)/Dt
                                end do
                            end if
                        end if
                        
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !  Apply moment
                        do K = 1,3
                            FM(K,I) = - rolling_moment(K) - twisting_moment(K) + FM(K,I) 
                        end do
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        
                        !  cohesive force
                        do K = 1,3
                            cohesive_force(K) = - m_c*(m_Beta*Rij)**2*DistU(K)
                        end do
                        do K = 1,3
                            F(K,I) = F(K,I) - cohesive_force(K)
                        end do
                        !  memory the contact in the Hertz linklist.
                        if (associated(Temp%prev)) then
                            !  Temp is in center of linklist!!!
                            if (Temp%No .EQ. J) then
                                !  Have contacted.
                                do K = 1,3
                                    Temp%Hertz(K) = tangential_force(K)
                                    Temp%Mrot(K) = rolling_moment(K)
                                    Temp%Mtwist(K) = twisting_moment(K)
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
                                TempH = Nodelink(J,tangential_force,rolling_moment,twisting_moment,&
                                & touching,slipping,rolling,twisting,Temp,Temp%next)
                                if (associated(Temp%next)) Temp%next%prev => TempH
                                Temp%next => TempH
                                Head(I)%No = LenNode + 1
                                Tail => TempH
                            end if
                        else
                            !  Temp is Head of linklist!!!
                            allocate(TempH)
                            TempH = Nodelink(J,tangential_force,rolling_moment,twisting_moment,&
                            & touching,slipping,rolling,twisting,Temp,Temp%next)
                            if (associated(Temp%next)) Temp%next%prev => TempH
                            Temp%next => TempH
                            Head(I)%No = LenNode + 1
                            Tail => TempH
                        end if
                    else
                        !  memory the separation in the Hertz linklist.
                        if (associated(Temp%prev) .AND. Temp%No.EQ.J) then
                            !  Temp is center of linklist!!!
                            Temp%prev%next => Temp%next
                            if(associated(Temp%next)) Temp%next%prev => Temp%prev
                            Head(I)%No = LenNode - 1
                            Tail => Temp%prev
                            deallocate(Temp)
                            !  When else Temp is Head of linklist!!!
                        else
                            Tail => Temp
                        end if
                        Rij = R(I)*R(J)/(R(I)+R(J))
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
                                F(K,I) = F(K,I) - cohesive_force(K)
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
                            !  apply force
                            do K = 1,3
                                F(K,I) = F(K,I) - cohesive_force(K)
                            end do        
                        end if
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    end if
                    !if (I.EQ.PP .OR. J.EQ.PP) then
                    !    write(123,*)
                    !end if
                end if
            end do
        end if

        !################         Wall Y1          ################### 
        if (Tag3(I).EQ.1) then
            Tail => Head(I)
            do J = 1,N
                if (Tag4(J).EQ.1) then   !  .OR. I.EQ.N
                    !  Initialize state params
                    do K = 1,3
                        Dist(K) = X(K,J) -  X(K,I)
                        H(K) = 0.0D0
                        Mr(K) = 0.0D0
                        Mt(K) = 0.0D0
                    end do
                    Dist(2) = Dist(2) - LenBoxY
                    shearPBC = X(1,J) + LenBoxY*gamma*time
                    shearPBC = shearPBC - ((shearPBC+0.5D0*LenBoxX) - MODULO((shearPBC+0.5D0*LenBoxX),LenBoxX))
                    Dist(1) = shearPBC - X(1,I)
                    touching = .false.
                    slipping = .false.
                    rolling = .false.
                    twisting = .false.
                    !  Distance vector
                    DistS = Dist(1)*Dist(1) + Dist(2)*Dist(2) + Dist(3)*Dist(3)
                    DistL = sqrt(DistS)
                    Dn = R(I) + R(J) - DistL
                    !  lookup this contact(or not) in the Hertz list
                    LenNode = Head(I)%No
                    Temp => Tail
                    do while(associated(Temp%next))
                        Temp => Temp%next
                        if (Temp%No .EQ. J) then
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
                        else if (Temp%No .GT. J) then
                            Temp => Temp%prev
                            exit
                        end if
                    end do
                    !  When collision calculate the repulsive restoring spring force which is generated along the normal and tangential according to Hooke's law
                    if (Dn .GT. 0.0D0) then
                        DistR = 1.0D0/DistL
                        !  calculate the normal vector
                        do K=1,3
                            DistU(K) = Dist(K)*DistR
                        end do
                        Ap = (R(I)*R(I)-R(J)*R(J)+DistS)/2.0D0*DistR
                        An = DistL-Ap
#ifdef HertzMindlinVisco                    
                        !  calculate material constant
                        Rij = R(I)*R(J)/(R(I)+R(J))
                        Mij = Body(I)*Body(J)/(Body(I)+Body(J))
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
                        Rij = R(I)*R(J)/(R(I)+R(J))
                        Mij = Body(I)*Body(J)/(Body(I)+Body(J))
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
                            Vrel(K) = Xdot(K,J) - Xdot(K,I)
                        end do
                        Vrel(1) = Vrel(1) + LenBoxY*gamma
                        !  negative rotate relative velocity
                        Vrot(1) = (DistU(2)*W(3,J) - DistU(3)*W(2,J))*An + (DistU(2)*W(3,I) - DistU(3)*W(2,I))*Ap
                        Vrot(2) = (DistU(3)*W(1,J) - DistU(1)*W(3,J))*An + (DistU(3)*W(1,I) - DistU(1)*W(3,I))*Ap
                        Vrot(3) = (DistU(1)*W(2,J) - DistU(2)*W(1,J))*An + (DistU(1)*W(2,I) - DistU(2)*W(1,I))*Ap
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
                        Energy(I) = Energy(I) + 0.4D0*Kn*(Dn**2)

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

                        if (slipping) then  !  Have slipped
                            if (DsL .GT. 1.0e-8) then  !  Still slipping
                                do K = 1,3
                                    tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    tangential_force(K) = 0.0D0  !  Particle J
                                end do
                                slipping = .false.
                            end if
                        else
                            if (tangential_forceL .GT. normal_forceL*m_mu_s) then  !  Slipping
                                slipping = .true.
                                if (DsL .GT. 1.0e-14) then
                                    do K = 1,3
                                        tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL
                                    end do
                                else
                                    do K = 1,3
                                        tangential_force(K) = 0.0D0
                                    end do
                                end if
                            else
                                slipping = .false.  !  Sticking
                            end if
                        end if

                        touching = .true.
                        !  Apply force
                        do K = 1,3
                            F(K,I) = - normal_force(K) - tangential_force(K) + F(K,I)
                        end do
                        !  Apply moment                     
                        FM(1,I) = - Ap*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,I) 
                        FM(2,I) = - Ap*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,I) 
                        FM(3,I) = - Ap*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,I)

                        !  rolling
                        do K = 1,3
                            Dtheta(K) = (W(K,I)-W(K,J))*Dt
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

                        if (rolling) then  !  Have rolled
                            if (DthetaRL .GT. 1.0e-8) then  !  Still slipping
                                do K = 1,3
                                    rolling_moment(K) = 2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    rolling_moment(K) = 0.0D0  !  Particle J
                                end do
                                rolling = .false.
                            end if
                        else
                            if (rolling_momentL .GT. 2.1D0*0.25D0*m_Beta*Rij*normal_forceL) then  !  Rolling
                                rolling = .true.
                                if (DthetaRL .GT. 1.0e-14) then
                                    do K = 1,3
                                        rolling_moment(K) = 2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL
                                    end do
                                else
                                    do K = 1,3
                                        rolling_moment(K) = 0.0D0
                                    end do
                                end if
                            else
                                rolling = .false.  !  Sticking
                                do K = 1,3
                                    rolling_moment(K) = rolling_moment(K) - Cr*DthetaR(K)/Dt
                                end do
                            end if
                        end if                                  

                        !  twisting moment of Particle J
                        do K = 1,3
                            twisting_moment(K) = Kt*DthetaT(K) + Mt(K)
                        end do
                        twisting_momentL = sqrt(twisting_moment(1)*twisting_moment(1) + twisting_moment(2)*twisting_moment(2) + twisting_moment(3)*twisting_moment(3))                

                        if (twisting) then  !  Have twisted
                            if (DthetaTL .GT. 1.0e-8) then  !  Still slipping
                                do K = 1,3
                                    twisting_moment(K) = 0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    twisting_moment(K) = 0.0D0  !  Particle J
                                end do
                                twisting = .false.
                            end if
                        else
                            if (twisting_momentL .GT. 0.65D0*m_mu_s*m_Beta*Rij*normal_forceL) then  !  Rolling
                                twisting = .true.
                                if (DthetaTL .GT. 1.0e-14) then
                                    do K = 1,3
                                        twisting_moment(K) = 0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL
                                    end do
                                else
                                    do K = 1,3
                                        twisting_moment(K) = 0.0D0
                                    end do
                                end if
                            else
                                twisting = .false.  !  Sticking
                                do K = 1,3
                                    twisting_moment(K) = twisting_moment(K) - Ct*DthetaT(K)/Dt
                                end do
                            end if
                        end if                                  

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !  Apply moment
                        do K = 1,3
                            FM(K,I) = - rolling_moment(K) - twisting_moment(K) + FM(K,I) 
                        end do
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !  cohesive force
                        do K = 1,3
                            cohesive_force(K) = - m_c*(m_Beta*Rij)**2*DistU(K)
                        end do
                        do K = 1,3
                            F(K,I) = F(K,I) - cohesive_force(K)
                        end do
                        !  memory the contact in the Hertz linklist.
                        if (associated(Temp%prev)) then
                            !  Temp is in center of linklist!!!
                            if (Temp%No .EQ. J) then
                                !  Have contacted.
                                do K = 1,3
                                    Temp%Hertz(K) = tangential_force(K)
                                    Temp%Mrot(K) = rolling_moment(K)
                                    Temp%Mtwist(K) = twisting_moment(K)
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
                                TempH = Nodelink(J,tangential_force,rolling_moment,twisting_moment,&
                                & touching,slipping,rolling,twisting,Temp,Temp%next)
                                if (associated(Temp%next)) Temp%next%prev => TempH
                                Temp%next => TempH
                                Head(I)%No = LenNode + 1
                                Tail => TempH
                            end if
                        else
                            !  Temp is Head of linklist!!!
                            allocate(TempH)
                            TempH = Nodelink(J,tangential_force,rolling_moment,twisting_moment,&
                            & touching,slipping,rolling,twisting,Temp,Temp%next)
                            if (associated(Temp%next)) Temp%next%prev => TempH
                            Temp%next => TempH
                            Head(I)%No = LenNode + 1
                            Tail => TempH
                        end if
                    else
                        !  memory the separation in the Hertz linklist.
                        if (associated(Temp%prev) .AND. Temp%No.EQ.J) then
                            !  Temp is center of linklist!!!
                            Temp%prev%next => Temp%next
                            if(associated(Temp%next)) Temp%next%prev => Temp%prev
                            Head(I)%No = LenNode - 1
                            Tail => Temp%prev
                            deallocate(Temp)
                            !  When else Temp is Head of linklist!!!
                        else
                            Tail => Temp
                        end if
                        Rij = R(I)*R(J)/(R(I)+R(J))
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
                                F(K,I) = F(K,I) - cohesive_force(K)
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
                            !  apply force
                            do K = 1,3
                                F(K,I) = F(K,I) - cohesive_force(K)
                            end do        
                        end if
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    end if
                end if
            end do
        end if

        !################         Wall Y2          ################### 
        if (Tag4(I).EQ.1) then
            Tail => Head(I)
            do J = 1,N
                if (Tag3(J).EQ.1) then   !  .OR. I.EQ.N
                    !  Initialize state params
                    do K = 1,3
                        Dist(K) = X(K,J) -  X(K,I)
                        H(K) = 0.0D0
                        Mr(K) = 0.0D0
                        Mt(K) = 0.0D0
                    end do
                    Dist(2) = Dist(2) + LenBoxY
                    shearPBC = X(1,J) - LenBoxY*gamma*time
                    shearPBC = shearPBC - ((shearPBC+0.5D0*LenBoxX) - MODULO((shearPBC+0.5D0*LenBoxX),LenBoxX))
                    Dist(1) = shearPBC - X(1,I)
                    touching = .false.
                    slipping = .false.
                    rolling = .false.
                    twisting = .false.
                    !  Distance vector
                    DistS = Dist(1)*Dist(1) + Dist(2)*Dist(2) + Dist(3)*Dist(3)
                    DistL = sqrt(DistS)
                    Dn = R(I) + R(J) - DistL
                    !  lookup this contact(or not) in the Hertz list
                    LenNode = Head(I)%No
                    Temp => Tail
                    do while(associated(Temp%next))
                        Temp => Temp%next
                        if (Temp%No .EQ. J) then
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
                        else if (Temp%No .GT. J) then
                            Temp => Temp%prev
                            exit
                        end if
                    end do
                    !  When collision calculate the repulsive restoring spring force which is generated along the normal and tangential according to Hooke's law
                    if (Dn .GT. 0.0D0) then
                        DistR = 1.0D0/DistL
                        !  calculate the normal vector
                        do K=1,3
                            DistU(K) = Dist(K)*DistR
                        end do
                        Ap = (R(I)*R(I)-R(J)*R(J)+DistS)/2.0D0*DistR
                        An = DistL-Ap
#ifdef HertzMindlinVisco                    
                        !  calculate material constant
                        Rij = R(I)*R(J)/(R(I)+R(J))
                        Mij = Body(I)*Body(J)/(Body(I)+Body(J))
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
                        Rij = R(I)*R(J)/(R(I)+R(J))
                        Mij = Body(I)*Body(J)/(Body(I)+Body(J))
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
                            Vrel(K) = Xdot(K,J) - Xdot(K,I)
                        end do
                        Vrel(1) = Vrel(1) - LenBoxY*gamma
                        !  negative rotate relative velocity
                        Vrot(1) = (DistU(2)*W(3,J) - DistU(3)*W(2,J))*An + (DistU(2)*W(3,I) - DistU(3)*W(2,I))*Ap
                        Vrot(2) = (DistU(3)*W(1,J) - DistU(1)*W(3,J))*An + (DistU(3)*W(1,I) - DistU(1)*W(3,I))*Ap
                        Vrot(3) = (DistU(1)*W(2,J) - DistU(2)*W(1,J))*An + (DistU(1)*W(2,I) - DistU(2)*W(1,I))*Ap
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
                        Energy(I) = Energy(I) + 0.4D0*Kn*(Dn**2)

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

                        if (slipping) then  !  Have slipped
                            if (DsL .GT. 1.0e-8) then  !  Still slipping
                                do K = 1,3
                                    tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    tangential_force(K) = 0.0D0  !  Particle J
                                end do
                                slipping = .false.
                            end if
                        else
                            if (tangential_forceL .GT. normal_forceL*m_mu_s) then  !  Slipping
                                slipping = .true.
                                if (DsL .GT. 1.0e-14) then
                                    do K = 1,3
                                        tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL
                                    end do
                                else
                                    do K = 1,3
                                        tangential_force(K) = 0.0D0
                                    end do
                                end if
                            else
                                slipping = .false.  !  Sticking
                            end if
                        end if

                        touching = .true.
                        !  Apply force
                        do K = 1,3
                            F(K,I) = - normal_force(K) - tangential_force(K) + F(K,I)
                        end do
                        !  Apply moment                     
                        FM(1,I) = - Ap*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,I) 
                        FM(2,I) = - Ap*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,I) 
                        FM(3,I) = - Ap*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,I)

                        !  rolling
                        do K = 1,3
                            Dtheta(K) = (W(K,I)-W(K,J))*Dt
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

                        if (rolling) then  !  Have rolled
                            if (DthetaRL .GT. 1.0e-8) then  !  Still slipping
                                do K = 1,3
                                    rolling_moment(K) = 2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    rolling_moment(K) = 0.0D0  !  Particle J
                                end do
                                rolling = .false.
                            end if
                        else
                            if (rolling_momentL .GT. 2.1D0*0.25D0*m_Beta*Rij*normal_forceL) then  !  Rolling
                                rolling = .true.
                                if (DthetaRL .GT. 1.0e-14) then
                                    do K = 1,3
                                        rolling_moment(K) = 2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL
                                    end do
                                else
                                    do K = 1,3
                                        rolling_moment(K) = 0.0D0
                                    end do
                                end if
                            else
                                rolling = .false.  !  Sticking
                                do K = 1,3
                                    rolling_moment(K) = rolling_moment(K) - Cr*DthetaR(K)/Dt
                                end do
                            end if
                        end if                                  

                        !  twisting moment of Particle J
                        do K = 1,3
                            twisting_moment(K) = Kt*DthetaT(K) + Mt(K)
                        end do
                        twisting_momentL = sqrt(twisting_moment(1)*twisting_moment(1) + twisting_moment(2)*twisting_moment(2) + twisting_moment(3)*twisting_moment(3))                

                        if (twisting) then  !  Have twisted
                            if (DthetaTL .GT. 1.0e-8) then  !  Still slipping
                                do K = 1,3
                                    twisting_moment(K) = 0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL  !  Particle J
                                end do
                            else  !  Approach sticking
                                do K = 1,3
                                    twisting_moment(K) = 0.0D0  !  Particle J
                                end do
                                twisting = .false.
                            end if
                        else
                            if (twisting_momentL .GT. 0.65D0*m_mu_s*m_Beta*Rij*normal_forceL) then  !  Rolling
                                twisting = .true.
                                if (DthetaTL .GT. 1.0e-14) then
                                    do K = 1,3
                                        twisting_moment(K) = 0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL
                                    end do
                                else
                                    do K = 1,3
                                        twisting_moment(K) = 0.0D0
                                    end do
                                end if
                            else
                                twisting = .false.  !  Sticking
                                do K = 1,3
                                    twisting_moment(K) = twisting_moment(K) - Ct*DthetaT(K)/Dt
                                end do
                            end if
                        end if                                 

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !  Apply moment
                        do K = 1,3
                            FM(K,I) = - rolling_moment(K) - twisting_moment(K) + FM(K,I) 
                        end do
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !  cohesive force
                        do K = 1,3
                            cohesive_force(K) = - m_c*(m_Beta*Rij)**2*DistU(K)
                        end do
                        do K = 1,3
                            F(K,I) = F(K,I) - cohesive_force(K)
                        end do
                        !  memory the contact in the Hertz linklist.
                        if (associated(Temp%prev)) then
                            !  Temp is in center of linklist!!!
                            if (Temp%No .EQ. J) then
                                !  Have contacted.
                                do K = 1,3
                                    Temp%Hertz(K) = tangential_force(K)
                                    Temp%Mrot(K) = rolling_moment(K)
                                    Temp%Mtwist(K) = twisting_moment(K)
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
                                TempH = Nodelink(J,tangential_force,rolling_moment,twisting_moment,&
                                & touching,slipping,rolling,twisting,Temp,Temp%next)
                                if (associated(Temp%next)) Temp%next%prev => TempH
                                Temp%next => TempH
                                Head(I)%No = LenNode + 1
                                Tail => TempH
                            end if
                        else
                            !  Temp is Head of linklist!!!
                            allocate(TempH)
                            TempH = Nodelink(J,tangential_force,rolling_moment,twisting_moment,&
                            & touching,slipping,rolling,twisting,Temp,Temp%next)
                            if (associated(Temp%next)) Temp%next%prev => TempH
                            Temp%next => TempH
                            Head(I)%No = LenNode + 1
                            Tail => TempH
                        end if
                    else
                        !  memory the separation in the Hertz linklist.
                        if (associated(Temp%prev) .AND. Temp%No.EQ.J) then
                            !  Temp is center of linklist!!!
                            Temp%prev%next => Temp%next
                            if(associated(Temp%next)) Temp%next%prev => Temp%prev
                            Head(I)%No = LenNode - 1
                            Tail => Temp%prev
                            deallocate(Temp)
                            !  When else Temp is Head of linklist!!!
                        else
                            Tail => Temp
                        end if
                        Rij = R(I)*R(J)/(R(I)+R(J))
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
                                F(K,I) = F(K,I) - cohesive_force(K)
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
                            !  apply force
                            do K = 1,3
                                F(K,I) = F(K,I) - cohesive_force(K)
                            end do        
                        end if
                    end if
                end if
            end do
        end if      
    end do
    !$OMP END PARALLEL DO

    return
    end
