    !********************************************************************
    !     DEMBody 4.1
    !     ***********
    !
    !     Force for all Particles.
    !     -------------------------------
    !     @Using violent traverse to search contact.
    !     @This model is used to checking Lattice method.
    !
    !     @Using rolling friction similar to LIGGGHTS
    !     @Using Hertz-Mindlin contact model similar to Wada
    !     @Using Cohesion model similar to Scheeres
    !     @Using damping model applicable for ice ball
    !
    !********************************************************************
    subroutine forceTraverse()

    use global
    use omp_lib
    implicit none

    real(kind=8)  Dist(3),DistS,DistL,DistR,DistU(3)
    real(kind=8)  Vrel(3),Vrot(3),Vtot(3),ERR,Vnor(3),Vtan(3)
    real(kind=8)  normal_force(3),normal_forceL
    real(kind=8)  tangential_force(3),tangential_forceL
    real(kind=8)  rolling_moment(3),rolling_momentL
    real(kind=8)  twisting_moment(3),twisting_momentL
    real(kind=8)  cohesive_force(3)
    real(kind=8)  gravity_force(3)
    real(kind=8)  Ap,An
    real(kind=8)  Rij,Mij,Iij
    real(kind=8)  Kn,Cn,Ks,Cs,Kr,Cr,Kt,Ct,lnCOR
    real(kind=8)  Dn,Ds(3),DsL,Dtheta(3),DthetaL,DthetaR(3),DthetaRL,DthetaT(3),DthetaTL
    real(kind=8)  H(3),Mr(3),Mt(3)
    real(kind=8)  RV(3)
    logical :: slipping,rolling,twisting    !  State of friction of T
    logical :: touching                     !  State of touch of T
    integer :: I,J,K,L,LenNode              !  Linklist
    integer :: num1,num2,limit              !  Neighbor
    type(Nodelink),pointer :: Temp          !  Temporary pointer
    type(Nodelink),pointer :: TempH         !  Contact pointer
    real(kind=8) ShearPBC                   !  shearPBC
    real(kind=8) accMag                     !  Magtitude of contact acceleration                
    integer :: region                       !  Upper limit of Loops 
    logical :: check                        !  Check for meshgrid

    integer :: particleI,particleJ          !  ID of two interacting particle
    real(8) :: Mass,MassCenter(3)           !  Mass and MassCenter of Gravity Lattice
    
    real(8) :: ostart,oend

    character(30) :: FileNameHead
    character(30) :: FileNameForce
    
    F = 0.0D0
    FM = 0.0D0
    Energy = 0.0D0
    
    !write(FileNameHead,'(F10.5)') Time
    !FileNameHead = trim(FileNameHead)//'Head.txt'
    !open(123,FILE=FileNameHead)

    !ostart = omp_get_wtime()
    !  Loop over all bodies through the NodeTree.
    !$OMP PARALLEL DO &
    !$OMP& PRIVATE(check,shearPBC,Temp,TempH,LenNode,num1,num2,limit,&
    !$OMP& I,J,K,L,Dist,DistS,DistL,DistR,DistU,Vrel,Vrot,Vtot,ERR,Vnor,Vtan,&
    !$OMP& normal_force,normal_forceL,tangential_force,tangential_forceL,&
    !$OMP& rolling_moment,rolling_momentL,twisting_moment,twisting_momentL,cohesive_force,gravity_force,Ap,An,Rij,Mij,Iij,&
    !$OMP& Kn,Cn,Ks,Cs,Kr,Cr,Kt,Ct,lnCOR,Dn,Ds,DsL,Dtheta,DthetaL,DthetaR,DthetaRL,DthetaT,DthetaTL,H,Mr,Mt,RV,&
    !$OMP& slipping,rolling,twisting,touching,&
    !$OMP& particleI,particleJ) SCHEDULE(DYNAMIC)
    !  Loop over all particles
    do I = 1,N
        do J = 1,N
            if (I .NE. J) then
                particleI = I
                particleJ = J

                !if (particleI.EQ.PP .OR. particleJ.EQ.PP) then
                !    write(123,'(F15.5,2X,2I5,2X)',advance='no') Time,particleI,particleJ
                !end if

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
                !  lookup this contact(or not) in the Hertz list
                LenNode = Head(particleI)%No
                Temp => Head(particleI)
                if (LenNode .NE. 0) then
                    Temp => Head(particleI)%next
                    do L = 1,LenNode
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
                        else if (Temp%No.LT.particleJ .AND. associated(Temp%next)) then
                            Temp => Temp%next
                        else if (Temp%No .GT. particleJ) then
                            Temp => Temp%prev
                            exit
                        end if
                    end do
                end if
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
                    Energy(particleI) = Energy(particleI) + 0.4D0*Kn*(Dn**2)
                    !Energy(particleJ) = Energy(particleJ) + 0.4D0*Kn*(Dn**2)

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
                                Temp%Hertz(K) = tangential_force(K)
                                Temp%Mrot(K) = rolling_moment(K)
                                Temp%Mtwist(K) = twisting_moment(K)
                            end do
                            Temp%is_touching = touching
                            Temp%is_slipping = slipping
                            Temp%is_rolling = rolling
                            Temp%is_twisting = twisting
                            Temp%recordTime = Time+Dt
                        else
                            !  First contacted.
                            allocate(TempH)
                            TempH = Nodelink(particleJ,Time+Dt,tangential_force,rolling_moment,twisting_moment,&
                            & touching,slipping,rolling,twisting,Temp,Temp%next)
                            if (associated(Temp%next)) Temp%next%prev => TempH
                            Temp%next => TempH
                            Head(particleI)%No = LenNode + 1
                        end if
                    else
                        !  Temp is Head of linklist!!!
                        allocate(TempH)
                        TempH = Nodelink(particleJ,Time+Dt,tangential_force,rolling_moment,twisting_moment,&
                        & touching,slipping,rolling,twisting,Temp,Temp%next)
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
                !if (particleI.EQ.PP .OR. particleJ.EQ.PP) then
                !    write(123,*)
                !end if
            end if
        end do    
    end do                  
    ! $OMP END PARALLEL DO
    !oend = omp_get_wtime()
    !write(*,*) "force",(oend-ostart)
    
    !write(*,*) "Time: ",Time
    !write(*,"(A12,3F15.5)") "Particles: ",(F(K,PP),K=1,3)

    !ostart = omp_get_wtime()
    !  calculate mirrored force if using PBC & Shear PBC
    if (isPeriodic) then
        call forceMirror
    end if
    !oend = omp_get_wtime()
    !write(*,*) "Mirror",(oend-ostart)
    
    !ostart = omp_get_wtime()
    !  calculate force of contactable walls if using contactable walls
    if (isContactWall) then
        call forceContactWalls
    end if
    !oend = omp_get_wtime()
    !write(*,*) "Contact walls", (oend-ostart)
    
    !ostart = omp_get_wtime()
    !  calculate force of bonded walls if using bonded walls
    if (isBondedWall) then                     
        call forceBondedWalls 
    end if
    !oend = omp_get_wtime()
    !write(*,*) "Bonded walls",(oend-ostart)

    !ostart = omp_get_wtime()
    !  calculate force of trimesh walls if using trimesh walls
    if (isTriMeshWall) then                     
        call forceTriMeshWalls 
    end if
    !oend = omp_get_wtime()
    !write(*,*) "Tri meshes",(oend-ostart)
    !write(*,"(A12,3F15.5)") "TriMesh: ",(F(K,PP),K=1,3)
    
    !ostart = omp_get_wtime()
    !  calculate force of funnel walls if using funnel walls
    if (isFunnelWall) then
        call forceFunnelWalls
    end if
    !oend = omp_get_wtime()
    !write(*,*) "Funnel walls",(oend-ostart)
    
    !ostart = omp_get_wtime()
    !  calculate force of gravity body if using GravBody
    if (isGravBody) then
        call forceGravBody
    end if
    !oend = omp_get_wtime()
    !write(*,*) "Grav bodies",(oend-ostart)

    !ostart = omp_get_wtime()
    !$OMP PARALLEL DO PRIVATE(I,K,accMag)
    do I = 1,N
        do K = 1,3
            F(K,I) = F(K,I)/Body(I)
            FM(K,I) = FM(K,I)/Inertia(I)
        end do
        accMag = sqrt(F(1,I)**2 + F(2,I)**2+F(3,I)**2)
        if (accMag .GE. MAX_ACC) then
            write(*,*) 'exceed Max Acc',I,accMag
            do K = 1,3
                F(K,I) = F(K,I)/accMag*MAX_ACC
            end do
        end if
        F(3,I) = F(3,I) - G
    end do
    !$OMP END PARALLEL DO
    !oend = omp_get_wtime()
    !write(*,*) 'accelerate', (oend-ostart)    

    !ostart = omp_get_wtime()
    !  calculate Planet force if in Planet system
    if (isPlanet) then
        call Planet
    end if
    !oend = omp_get_wtime()
    !write(*,*) 'Planet', (oend-ostart)    
    
    !ostart = omp_get_wtime()
    !  calculate Noninertial force if in Rotary system
    if (isRotSystem) then
        call forceRotSystem
    end if
    !oend = omp_get_wtime()
    !write(*,*) 'Rot system', (oend-ostart)    
    
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
