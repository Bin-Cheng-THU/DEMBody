    !********************************************************************
    !     DEMBody 7.2
    !     ***********
    !
    !     Force for BiDisperse Gravity Body.
    !     --------------------------
    !      
    !     @Using rolling friction similar to LIGGGHTS
    !     @Using Hertz-Mindlin contact model similar to Wada
    !     @Using Cohesion model similar to Scheeres
    !     @Using damping model applicable for ice ball
    !
    !********************************************************************
    subroutine forceBiDisperseGravBody()

    use global
    use omp_lib
    implicit none

    real(kind=8)  Dist(3),DistS,DistL,DistR,DistU(3)
    real(kind=8)  Vrel(3),Vrot(3),Vtot(3),ERR,Vnor(3),Vtan(3)
    real(kind=8)  normal_force(3),normal_forceL
    real(kind=8)  tangential_force(3),tangential_forceL,tangential_history(3)
    real(kind=8)  rolling_moment(3),rolling_momentL,rolling_history(3)
    real(kind=8)  twisting_moment(3),twisting_momentL,twisting_history(3)
    real(kind=8)  cohesive_force(3)
    real(kind=8)  gravity_force(3)
    real(kind=8)  Ap,An
    real(kind=8)  Rij,Mij,Iij
    real(kind=8)  Kn,Cn,Ks,Cs,Kr,Cr,Kt,Ct,lnCOR
    real(kind=8)  Dn,Ds(3),DsL,Dtheta(3),DthetaL,DthetaR(3),DthetaRL,DthetaT(3),DthetaTL
    real(kind=8)  H(3),Mr(3),Mt(3)
    real(kind=8)  RV(3)
    logical :: slipping,rolling,twisting             !  State of friction of T
    logical :: touching                              !  State of touch of T
    integer :: I,J,K,L,II,LenNode                    !  Iterator
    type(Nodelink),pointer :: Temp                   !  Temporary pointer
    type(Nodelink),pointer :: TempH                  !  Contact pointer
    real(kind=8)  rPanV(3),rPan,center,rOrigM

    real(8) :: ostart,oend

    gravBodyF = 0.0D0
    gravBodyFM = 0.0D0
        
    !ostart = omp_get_wtime()
    !  Loop over all bodies.
    !$OMP PARALLEL DO REDUCTION(+:gravBodyF) REDUCTION(+:gravBodyFM)&
    !$OMP& PRIVATE(Temp,TempH,LenNode,&
    !$OMP& I,J,K,L,II,Dist,DistS,DistL,DistR,DistU,Vrel,Vrot,Vtot,ERR,Vnor,Vtan,&
    !$OMP& normal_force,normal_forceL,tangential_force,tangential_forceL,tangential_history,&
    !$OMP& rolling_moment,rolling_momentL,rolling_history,twisting_moment,twisting_momentL,twisting_history,cohesive_force,gravity_force,Ap,An,Rij,Mij,Iij,&
    !$OMP& Kn,Cn,Ks,Cs,Kr,Cr,Kt,Ct,lnCOR,Dn,Ds,DsL,Dtheta,DthetaL,DthetaR,DthetaRL,DthetaT,DthetaTL,H,Mr,Mt,RV,&
    !$OMP& slipping,rolling,twisting,touching) SCHEDULE(GUIDED)
    do I = 1,biDisperseNum
        do K = 1,3
            Dist(K) = gravBodyX(K) - biDisperseX(K,I)
            H(K) = 0.0D0
            Mr(K) = 0.0D0
            Mt(K) = 0.0D0
        end do
        
        !  Distance vector
        DistS = Dist(1)*Dist(1) + Dist(2)*Dist(2) + Dist(3)*Dist(3)
        DistL = sqrt(DistS)
        DistR = 1.0D0/DistL

        do K = 1,3
            gravity_force(K) = GravConst*biDisperseBody(I)*gravBodyBody*Dist(K)*DistR*DistR*DistR
        end do

        do K = 1,3
            biDisperseF(K,I) = biDisperseF(K,I) + gravity_force(K)
            gravBodyF(K) = gravBodyF(K) - gravity_force(K)
        end do

        ERR = DistL - biDisperseR(I) - gravBodyR
        if (ERR .LE. 2.0D0*biDisperseScale*LatDx) then
            touching = .false.
            slipping = .false.
            rolling = .false.        
            twisting = .false.
            
            Dn = biDisperseR(I) + gravBodyR - DistL
            !  lookup this contact(or not) in the Hertz list
            LenNode = HeadBiDisperse(I)%No
            Temp => HeadBiDisperse(I)
            if (LenNode .NE. 0) then
                Temp => HeadBiDisperse(I)%next
                do L = 1,LenNode
                    if (Temp%No .EQ. gravBodyTag+biDisperseNum) then
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
                    else if (Temp%No.LT.gravBodyTag+biDisperseNum .AND. associated(Temp%next)) then
                        Temp => Temp%next
                    else if (Temp%No .GT. gravBodyTag+biDisperseNum) then
                        Temp => Temp%prev
                        exit
                    end if
                end do
            end if
            !  When collision calculate the repulsive restoring spring force which is generated along the normal and tangential according to Hertzian law
            if (Dn .GT. 0.0D0) then 
                !  calculate the normal vector
                do K=1,3
                    DistU(K) = Dist(K)*DistR
                end do
                Ap = (biDisperseR(I)*biDisperseR(I)-gravBodyR*gravBodyR+DistS)/2.0D0*DistR
                An = DistL-Ap
#ifdef HertzMindlinVisco
                !  calculate material constant
                Rij = biDisperseR(I)*gravBodyR/(biDisperseR(I)+gravBodyR)
                Mij = biDisperseBody(I)*gravBodyBody/(biDisperseBody(I)+gravBodyBody)
                Kn = 2.0D0*m_Ei*sqrt(Rij*Dn)/(3.0D0*(1.0D0-m_nu*m_nu))
                Cn = -Kn*m_A
                Ks = 2.0D0*m_Ei/(1.0D0+m_nu)/(2.0D0-m_nu)*sqrt(Rij)*sqrt(Dn)
                !  select tangential damping mode
                if (m_COR > 1.0D0) then
                    Cs = -2.0D0*m_Ei/(1.0D0+m_nu)/(2.0D0-m_nu)*sqrt(Dn)*m_A
                elseif (m_COR >= 0.0D0) then
                    lnCOR=log(m_COR)
                    Cs = 2.0D0*sqrt(5.0D0/6.0D0)*lnCOR/sqrt(lnCOR**2+3.1415926D0**2) &
                    & *sqrt(2.0D0*Mij*m_Ei/(1.0D0+m_nu)/(2.0D0-m_nu))*(Rij**0.25)*(Dn**0.25)
                else
                    Cs = 0.0D0
                end if
                Kr = 0.25D0*Kn*(m_Beta*Rij)**2
                Cr = 0.25D0*Cn*(m_Beta*Rij)**2
                Kt = 0.5D0*Ks*(m_Beta*Rij)**2
                Ct = 0.5D0*Cs*(m_Beta*Rij)**2
#elif HertzMindlinResti
                !  calculate material constant
                Rij = biDisperseR(I)*gravBodyR/(biDisperseR(I)+gravBodyR)
                Mij = biDisperseBody(I)*gravBodyBody/(biDisperseBody(I)+gravBodyBody)
                Kn = 2.0D0*m_Ei*sqrt(Rij*Dn)/(3.0D0*(1.0D0-m_nu*m_nu))
                lnCOR = log(m_COR)
                Cn = 2.0D0*sqrt(5.0D0/6.0D0)*lnCOR/sqrt(lnCOR**2+3.1415926D0**2) &
                & *sqrt(Mij*m_Ei/(1.0D0-m_nu*m_nu))*(Rij**0.25)*(Dn**0.25)
                Ks = 2.0D0*m_Ei/(1.0D0+m_nu)/(2.0D0-m_nu)*sqrt(Rij)*sqrt(Dn)
                Cs = 2.0D0*sqrt(5.0D0/6.0D0)*lnCOR/sqrt(lnCOR**2+3.1415926D0**2) &
                & *sqrt(2.0D0*Mij*m_Ei/(1.0D0+m_nu)/(2.0D0-m_nu))*(Rij**0.25)*(Dn**0.25)
                Kr = 0.25D0*Kn*(m_Beta*Rij)**2
                Cr = 0.25D0*Cn*(m_Beta*Rij)**2
                Kt = 0.5D0*Ks*(m_Beta*Rij)**2
                Ct = 0.5D0*Cs*(m_Beta*Rij)**2
#endif
                !  translate relative velocity
                do K = 1,3
                    Vrel(K) = gravBodyXdot(K) - biDisperseXdot(K,I)
                end do
                !  rotate relative velocity
                Vrot(1) = (DistU(2)*gravBodyW(3) - DistU(3)*gravBodyW(2))*An + (DistU(2)*biDisperseW(3,I) - DistU(3)*biDisperseW(2,I))*Ap
                Vrot(2) = (DistU(3)*gravBodyW(1) - DistU(1)*gravBodyW(3))*An + (DistU(3)*biDisperseW(1,I) - DistU(1)*biDisperseW(3,I))*Ap
                Vrot(3) = (DistU(1)*gravBodyW(2) - DistU(2)*gravBodyW(1))*An + (DistU(1)*biDisperseW(2,I) - DistU(2)*biDisperseW(1,I))*Ap
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

                !  normal force of Gravity Body
                do K = 1,3
                    normal_force(K) = Kn*Dn*DistU(K) + Cn*Vnor(K)
                end do
                normal_forceL = sqrt(normal_force(1)*normal_force(1) + normal_force(2)*normal_force(2) + normal_force(3)*normal_force(3))

                !  Add energy
                biDisperseEnergy(I) = biDisperseEnergy(I) + 0.2D0*Kn*(Dn**2)

                !  tangential deform
                do K = 1,3
                    Ds(K) = Vtan(K)*Dt
                end do
                DsL = sqrt(Ds(1)*Ds(1) + Ds(2)*Ds(2) + Ds(3)*Ds(3))

                !  tangential force of Gravity Body
                do K = 1,3
                    tangential_force(K) = - Ks*Ds(K) + Cs*Vtan(K) + H(K)
                end do
                tangential_forceL = sqrt(tangential_force(1)*tangential_force(1) + tangential_force(2)*tangential_force(2) + tangential_force(3)*tangential_force(3))

                !if (slipping) then  !  Have slipped
                !    if (DsL .GT. 1.0e-14) then  !  Still slipping
                !        do K = 1,3
                !            tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL  !  Gravity Body
                !        end do
                !    else  !  Approach sticking
                !        do K = 1,3
                !            tangential_force(K) = 0.0D0  !  Gravity Body
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
                            biDisperseEnergy(I) = biDisperseEnergy(I) + 0.5D0*m_mu_d*normal_forceL*DsL
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
                    gravBodyF(K) = normal_force(K) + tangential_force(K) + gravBodyF(K)
                    biDisperseF(K,I) = - normal_force(K) - tangential_force(K) + biDisperseF(K,I)
                end do
                !  Apply moment
                gravBodyFM(1) = - An*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + gravBodyFM(1) 
                gravBodyFM(2) = - An*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + gravBodyFM(2) 
                gravBodyFM(3) = - An*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + gravBodyFM(3)                     
                biDisperseFM(1,I) = - Ap*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + biDisperseFM(1,I) 
                biDisperseFM(2,I) = - Ap*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + biDisperseFM(2,I) 
                biDisperseFM(3,I) = - Ap*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + biDisperseFM(3,I)

                !  rolling
                do K = 1,3
                    Dtheta(K) = (biDisperseW(K,I)-gravBodyW(K))*Dt
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

                !  rolling moment of Gravity Body
                do K = 1,3
                    rolling_moment(K) = Kr*DthetaR(K) + Mr(K)
                end do
                rolling_momentL = sqrt(rolling_moment(1)*rolling_moment(1) + rolling_moment(2)*rolling_moment(2) + rolling_moment(3)*rolling_moment(3))                

                !if (rolling) then  !  Have rolled
                !    if (DthetaRL .GT. 1.0e-14) then  !  Still slipping
                !        do K = 1,3
                !            rolling_moment(K) = 2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL  !  Gravity Body
                !        end do
                !    else  !  Approach sticking
                !        do K = 1,3
                !            rolling_moment(K) = 0.0D0  !  Gravity Body
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
                        else
                            do K = 1,3
                                rolling_moment(K) = 0.0D0
                                rolling_history(K) = rolling_moment(K)
                            end do
                        end if
                    else
                        rolling = .false. !  Sticking
                        do K = 1,3
                            rolling_history(K) = rolling_moment(K)
                            rolling_moment(K) = rolling_moment(K) - Cr*DthetaR(K)/Dt
                        end do
                    end if
                !end if                                  

                !  twisting moment of Gravity Body
                do K = 1,3
                    twisting_moment(K) = Kt*DthetaT(K) + Mt(K)
                end do
                twisting_momentL = sqrt(twisting_moment(1)*twisting_moment(1) + twisting_moment(2)*twisting_moment(2) + twisting_moment(3)*twisting_moment(3))                

                !if (twisting) then  !  Have twisted
                !    if (DthetaTL .GT. 1.0e-14) then  !  Still slipping
                !        do K = 1,3
                !            twisting_moment(K) = 0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL  !  Gravity Body
                !        end do
                !    else  !  Approach sticking
                !        do K = 1,3
                !            twisting_moment(K) = 0.0D0  !  Gravity Body
                !        end do
                !        twisting = .false.
                !    end if
                !else
                    if (twisting_momentL .GT. 0.65D0*m_mu_s*m_Beta*Rij*normal_forceL) then  !  Rolling
                        twisting = .true.
                        if (DthetaTL .GT. 1.0e-14) then
                            do K = 1,3
                                twisting_moment(K) = 0.65D0*m_mu_s*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL
                                twisting_history(K) = twisting_moment(K)
                            end do
                        else
                            do K = 1,3
                                twisting_moment(K) = 0.0D0
                                twisting_history(K) = twisting_moment(K)
                            end do
                        end if
                    else
                        twisting = .false. !  Sticking
                        do K = 1,3
                            twisting_history(K) = twisting_moment(K)
                            twisting_moment(K) = twisting_moment(K) - Cr*DthetaT(K)/Dt
                        end do
                    end if
                !end if   
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
                !  Apply moment
                do K = 1,3
                    gravBodyFM(K) = + rolling_moment(K) + twisting_moment(K) + gravBodyFM(K) 
                    biDisperseFM(K,I) = - rolling_moment(K) - twisting_moment(K) + biDisperseFM(K,I) 
                end do
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

                !  cohesive force
                do K = 1,3
                    cohesive_force(K) = - m_c*(m_Beta*Rij)**2*DistU(K)
                end do
                do K = 1,3
                    gravBodyF(K) = gravBodyF(K) + cohesive_force(K)
                    biDisperseF(K,I) = biDisperseF(K,I) - cohesive_force(K)
                end do

                !  memory the contact in the Hertz linklist.
                if (associated(Temp%prev)) then
                    !  Temp is in center of linklist!!!
                    if (Temp%No .EQ. gravBodyTag+biDisperseNum) then
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
                    else
                        !  First contacted.
                        allocate(TempH)
                        TempH = Nodelink(gravBodyTag+biDisperseNum,tangential_history,rolling_history,twisting_history,&
                        & touching,slipping,rolling,twisting,Temp,Temp%next)
                        if (associated(Temp%next)) Temp%next%prev => TempH
                        Temp%next => TempH
                        HeadBiDisperse(I)%No = LenNode + 1
                    end if
                else
                    !  Temp is Head of linklist!!!
                    allocate(TempH)
                    TempH = Nodelink(gravBodyTag+biDisperseNum,tangential_history,rolling_history,twisting_history,&
                    & touching,slipping,rolling,twisting,Temp,Temp%next)
                    if (associated(Temp%next)) Temp%next%prev => TempH
                    Temp%next => TempH
                    HeadBiDisperse(I)%No = LenNode + 1
                end if
            else
                !  memory the separation in the Hertz linklist.
                if (associated(Temp%prev) .AND. Temp%No.EQ.gravBodyTag+biDisperseNum) then
                    !  Temp is center of linklist!!!
                    Temp%prev%next => Temp%next
                    if(associated(Temp%next)) Temp%next%prev => Temp%prev
                    HeadBiDisperse(I)%No = LenNode - 1
                    deallocate(Temp)
                    !  When else Temp is Head of linklist!!!
                end if
                Rij = biDisperseR(I)*gravBodyR/(biDisperseR(I)+gravBodyR)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (Dn > -Rij*(m_r_cut-1.0D0)*0.5D0) then
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
                        gravBodyF(K) = gravBodyF(K) + cohesive_force(K)
                        biDisperseF(K,I) = biDisperseF(K,I) - cohesive_force(K)
                    end do                         
                else if (Dn > -Rij*(m_r_cut-1.0D0)) then
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
                        gravBodyF(K) = gravBodyF(K) + cohesive_force(K)
                        biDisperseF(K,I) = biDisperseF(K,I) - cohesive_force(K)
                    end do        
                end if                
            end if
        end if
    end do
    !$OMP END PARALLEL DO
    !oend = omp_get_wtime()
    !write(*,*) "BiDisperse GravBody",(oend-ostart)
    
    !*********modification for Phobos
    do K = 1,3
        gravBodyF(K) = 0.0!gravBodyF(K)/gravBodyBody
        gravBodyFM(K) = 0.0!gravBodyFM(K)/gravBodyInertia
    end do
        
    return
    end
