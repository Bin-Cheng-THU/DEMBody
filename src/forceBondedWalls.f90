    !********************************************************************
    !     DEMBody 3.0
    !     ***********
    !
    !     Force for bonded walls.
    !     --------------------------
    !      
    !     @Using rolling friction similar to LIGGGHTS
    !     @Using Hertz-Mindlin contact model similar to Wada
    !     @No Cohesion model in walls
    !     @No Rolling model in bonded walls
    !     @Using damping model applicable for ice ball
    !
    !********************************************************************
    subroutine forceBondedWalls()

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
    logical :: slipping,rolling       !  State of friction of T
    logical :: touching               !  State of touch of T
    integer :: I,J,K,L,II,LenNode     !  Iterator
    type(Nodelink),pointer :: Temp    !  Temporary pointer
    type(Nodelink),pointer :: TempH   !  Contact pointer
    integer :: OMP_bondedWallTag      !  Tag for walls in OMP
    real(kind=8) OMP_bondedWallPoint(3),OMP_bondedWallVectorN(3)  !  Point and Normal Vector for walls in OMP
    real(kind=8) OMP_bondedWallVectorTx(3),OMP_bondedWallVectorTy(3)  !  Tangential Vector for walls in OMP
    real(kind=8) OMP_bondedWallLx,OMP_bondedWallLy  !  Length for walls in OMP
    logical :: enterFlag
    real(kind=8)  RVb(3),RVx,RVy,RVn(3),RVt(3),norm
    real(kind=8)  RVo(3),Vwall(3)
    real(kind=8)  temp_bondedWallFM(3)

    bondedWallF = 0.0D0
    bondedWallFM = 0.0D0
    
    do II = 1,bondedWallNum
        !  Allocate to variables in OMP (stack memory) and refresh Points & Vectors of Bonded Walls
        OMP_bondedWallTag = bondedWallTag(II)
        do K = 1,3
            OMP_bondedWallPoint(K) = bondedWallX(K) + bondedWallMatB(K,1)*bondedWallPoint(1,II) + bondedWallMatB(K,2)*bondedWallPoint(2,II) + bondedWallMatB(K,3)*bondedWallPoint(3,II)
            OMP_bondedWallVectorN(K) = bondedWallMatB(K,1)*bondedWallVectorN(1,II) + bondedWallMatB(K,2)*bondedWallVectorN(2,II) + bondedWallMatB(K,3)*bondedWallVectorN(3,II)
            OMP_bondedWallVectorTx(K) = bondedWallMatB(K,1)*bondedWallVectorTx(1,II) + bondedWallMatB(K,2)*bondedWallVectorTx(2,II) + bondedWallMatB(K,3)*bondedWallVectorTx(3,II)
            OMP_bondedWallVectorTy(K) = bondedWallMatB(K,1)*bondedWallVectorTy(1,II) + bondedWallMatB(K,2)*bondedWallVectorTy(2,II) + bondedWallMatB(K,3)*bondedWallVectorTy(3,II)
        end do
        OMP_bondedWallLx = bondedWallLx(II)
        OMP_bondedWallLy = bondedWallLy(II)
        
        !  Loop over all bodies.
        !$OMP PARALLEL DO REDUCTION(+:F) REDUCTION(+:FM) &
        !$OMP& REDUCTION(+:Energy) &
        !$OMP& REDUCTION(+:bondedWallF) REDUCTION(+:bondedWallFM) &
        !$OMP& firstprivate(OMP_bondedWallTag,OMP_bondedWallPoint,OMP_bondedWallVectorN,OMP_bondedWallVectorTx,OMP_bondedWallVectorTy) &
        !$OMP& PRIVATE(RVb,RVx,RVy,RVn,RVt,RVo,Vwall,norm,enterFlag)&
        !$OMP& PRIVATE(Temp,TempH,LenNode,&
        !$OMP& I,J,K,L,II,Dist,DistS,DistL,DistR,DistU,Vrel,Vrot,Vtot,ERR,Vnor,Vtan,&
        !$OMP& normal_force,normal_forceL,tangential_force,tangential_forceL,&
        !$OMP& rolling_moment,rolling_momentL,cohesive_force,Ap,An,Rij,Mij,Iij,&
        !$OMP& Kn,Cn,Kt,Ct,Kr,Cr,lnCOR,Dn,Ds,DsL,Dtheta,DthetaL,H,M,RV,&
        !$OMP& slipping,rolling,touching) SCHEDULE(DYNAMIC)
        do I = 1,N
            enterFlag = .False.
            do K = 1,3
                RV(K) = X(K,I) - OMP_bondedWallPoint(K)
            end do
            ERR = RV(1)*OMP_bondedWallVectorN(1) + RV(2)*OMP_bondedWallVectorN(2) + RV(3)*OMP_bondedWallVectorN(3)
            if (ERR > 0.0 .AND. ERR < Dx) then
                do K = 1,3
                    RVb(K) = RV(K) - ERR*OMP_bondedWallVectorN(K)
                end do 
                RVx = RVb(1)*OMP_bondedWallVectorTx(1) + RVb(2)*OMP_bondedWallVectorTx(2) + RVb(3)*OMP_bondedWallVectorTx(3)
                RVy = RVb(1)*OMP_bondedWallVectorTy(1) + RVb(2)*OMP_bondedWallVectorTy(2) + RVb(3)*OMP_bondedWallVectorTy(3)
                if (ABS(RVx).LE.(0.5*OMP_bondedWallLx) .AND. ABS(RVy).LE.(0.5*OMP_bondedWallLy)) then
                    enterFlag = .True.
                    do K = 1,3
                        RVt(K) = RVb(K)
                        RVn(K) = ERR*OMP_bondedWallVectorN(K)
                    end do
                else if (ABS(RVx).LE.(0.5*OMP_bondedWallLx) .AND. ABS(RVy).LE.(0.5*OMP_bondedWallLy+Dx)) then
                    enterFlag = .True.
                    do K = 1,3
                        RVt(K) = RVy/ABS(RVy)*0.5*OMP_bondedWallLy*OMP_bondedWallVectorTy(K) + RVx*OMP_bondedWallVectorTx(K)
                        RVn(K) = RV(K) - RVt(K)
                    end do
                else if (ABS(RVx).LE.(0.5*OMP_bondedWallLx+Dx) .AND. ABS(RVy).LE.(0.5*OMP_bondedWallLy)) then
                    enterFlag = .True.
                    do K = 1,3
                        RVt(K) = RVx/ABS(RVx)*0.5*OMP_bondedWallLx*OMP_bondedWallVectorTx(K) + RVy*OMP_bondedWallVectorTy(K)
                        RVn(K) = RV(K) - RVt(K)
                    end do                        
                else if (ABS(RVx).LE.(0.5*OMP_bondedWallLx+Dx) .AND. ABS(RVy).LE.(0.5*OMP_bondedWallLy+Dx)) then
                    enterFlag = .True.
                    do K = 1,3
                        RVt(K) = 0.5*OMP_bondedWallLx*OMP_bondedWallVectorTx(K)*RVx/ABS(RVx) + 0.5*OMP_bondedWallLy*OMP_bondedWallVectorTy(k)*RVy/ABS(RVy)
                        RVn(K) = RV(K) - RVt(K)
                    end do     
                end if
            end if
            
            if (enterFlag) then      
                !  normal vector
                do K = 1,3
                    Dist(K) = RVn(K)
                    H(K) = 0.0
                    M(K) = 0.0
                end do
                touching = .false.
                slipping = .false.
                rolling = .false.            
                !  Distance vector
                DistS = Dist(1)*Dist(1) + Dist(2)*Dist(2) + Dist(3)*Dist(3)
                DistL = sqrt(DistS)
                Dn = R(I) - DistL
                !  lookup this contact(or not) in the Hertz list
                LenNode = Head(I)%No
                Temp => Head(I)
                if (LenNode .NE. 0) then
                    Temp => Head(I)%next
                    do L = 1,LenNode
                        if (Temp%No .EQ. OMP_bondedWallTag) then
                            do K = 1,3
                                H(K) = Temp%Hertz(K)
                                M(K) = Temp%Mrot(K)
                            end do
                            touching = Temp%is_touching
                            slipping = Temp%is_slipping
                            rolling = Temp%is_rolling
                            exit
                        else if (Temp%No.LT.OMP_bondedWallTag .AND. associated(Temp%next)) then
                            Temp => Temp%next
                        else if (Temp%No .GT. OMP_bondedWallTag) then
                            Temp => Temp%prev
                            exit
                        end if
                    end do
                end if
                !  When collision calculate the repulsive restoring spring force which is generated along the normal and tangential according to Hertzian law
                if (Dn .GT. 0.0) then 

                    DistR = 1.0/DistL
                    !  calculate the normal vector
                    do K=1,3
                        DistU(K) = Dist(K)*DistR
                    end do
                    Ap = DistL
#ifdef HertzMindlinVisco                    
                    !  calculate material constant
                    Rij = R(I)
                    Mij = Body(I)
                    Iij = 3.5*Inertia(I)
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
                    Rij = R(I)
                    Mij = Body(I)
                    Iij = 3.5*Inertia(I)
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
                    !  distance vector relative to Bonded Wall Center
                    do K = 1,3
                        RVo(K) = RVt(K) + OMP_bondedWallPoint(K) - bondedWallX(K)
                    end do
                    !  rotate velocity of contact point
                    Vwall(1) = bondedWallW(2)*RVo(3) - bondedWallW(3)*RVo(2)
                    Vwall(2) = bondedWallW(3)*RVo(1) - bondedWallW(1)*RVo(3)
                    Vwall(3) = bondedWallW(1)*RVo(2) - bondedWallW(2)*RVo(1)                    
                    !  translate relative velocity
                    do K = 1,3
                        Vrel(K) = Xdot(K,I) - bondedWallXdot(K)
                    end do
                    !  rotate relative velocity
                    Vrot(1) = (DistU(2)*W(3,I) - DistU(3)*W(2,I))*Ap - Vwall(1)
                    Vrot(2) = (DistU(3)*W(1,I) - DistU(1)*W(3,I))*Ap - Vwall(2)
                    Vrot(3) = (DistU(1)*W(2,I) - DistU(2)*W(1,I))*Ap - Vwall(3)
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

                    !  normal force of Particle I
                    do K = 1,3
                        normal_force(K) = Kn*(Dn**1.5)*DistU(K) + Cn*Vnor(K)
                    end do
                    normal_forceL = sqrt(normal_force(1)*normal_force(1) + normal_force(2)*normal_force(2) + normal_force(3)*normal_force(3))

                    !  Add energy
                    Energy(I) = Energy(I) + 0.2*Kn*(Dn**2.5)
                    
                    !  tangential deform
                    do K = 1,3
                        Ds(K) = Vtan(K)*Dt
                    end do
                    DsL = sqrt(Ds(1)*Ds(1) + Ds(2)*Ds(2) + Ds(3)*Ds(3))

                    !  tangential force of Particle I
                    do K = 1,3
                        tangential_force(K) = - Kt*Ds(K) + Ct*Vtan(K) + H(K)
                    end do
                    tangential_forceL = sqrt(tangential_force(1)*tangential_force(1) + tangential_force(2)*tangential_force(2) + tangential_force(3)*tangential_force(3))

                    if (slipping) then  !  Have slipped
                        if (DsL .GT. 1e-8) then  !  Still slipping
                            do K = 1,3
                                tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL  !  Particle I
                            end do
                        else  !  Approach sticking
                            do K = 1,3
                                tangential_force(K) = 0.0  !  Particle I
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
                        F(K,I) = normal_force(K) + tangential_force(K) + F(K,I)
                        bondedWallF(K) = -normal_force(K) - tangential_force(K) + bondedWallF(K)
                    end do
                    !  Apply moment                  
                    FM(1,I) = - Ap*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,I) 
                    FM(2,I) = - Ap*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,I) 
                    FM(3,I) = - Ap*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,I) 
                    bondedWallFM(1) = - RVo(2)*(normal_force(3) + tangential_force(3)) + RVo(3)*(normal_force(2) + tangential_force(2)) + bondedWallFM(1)
                    bondedWallFM(2) = - RVo(3)*(normal_force(1) + tangential_force(1)) + RVo(1)*(normal_force(3) + tangential_force(3)) + bondedWallFM(2)
                    bondedWallFM(3) = - RVo(1)*(normal_force(2) + tangential_force(2)) + RVo(2)*(normal_force(1) + tangential_force(1)) + bondedWallFM(3)       
                    
                    do K = 1,3
                        rolling_moment(K) = 0.0
                    end do
                    
                    !!  rolling deform
                    !do K = 1,3
                    !    Dtheta(K) = (W(K,I))*Dt
                    !end do
                    !DthetaL = sqrt(Dtheta(1)*Dtheta(1) + Dtheta(2)*Dtheta(2) + Dtheta(3)*Dtheta(3))
                    !
                    !!  rolling moment of Particle I
                    !do K = 1,3
                    !    rolling_moment(K) = - Kr*Dtheta(K) + M(K)
                    !end do
                    !rolling_momentL = sqrt(rolling_moment(1)*rolling_moment(1) + rolling_moment(2)*rolling_moment(2) + rolling_moment(3)*rolling_moment(3))                
                    !
                    !if (rolling) then  !  Have rolled
                    !    if (DthetaL .GT. 1e-8) then  !  Still slipping
                    !        do K = 1,3
                    !            rolling_moment(K) = -m_mu_r*Rij*normal_forceL*Dtheta(K)/DthetaL  !  Particle I
                    !        end do
                    !    else  !  Approach sticking
                    !        do K = 1,3
                    !            rolling_moment(K) = 0.0  !  Particle I
                    !        end do
                    !        rolling = .false.
                    !    end if
                    !else
                    !    if (rolling_momentL .GT. normal_forceL*Rij*m_mu_r) then  !  Rolling
                    !        rolling = .true.
                    !        if (DthetaL .GT. 1e-14) then
                    !            do K = 1,3
                    !                rolling_moment(K) = -m_mu_r*Rij*normal_forceL*Dtheta(K)/DthetaL
                    !            end do
                    !        else
                    !            do K = 1,3
                    !                rolling_moment(K) = 0.0
                    !            end do
                    !        end if
                    !    else
                    !        rolling = .false. !  Sticking
                    !        do K = 1,3
                    !            rolling_moment(K) = rolling_moment(K) - Cr*(W(K,I))
                    !        end do
                    !    end if
                    !end if                                  
                    !
                    !!  Apply moment
                    !do K = 1,3
                    !    FM(K,I) = rolling_moment(K) + FM(K,I)
                    !end do

                    !!  cohesive force
                    !do K = 1,3
                    !    cohesive_force(K) = - 1.0*m_Beta*3.1415926*Rij**2*DistU(K)
                    !end do
                    !do K = 1,3
                    !    F(K,I) = F(K,I) + cohesive_force(K)
                    !end do
                    
                    !  memory the contact in the Hertz linklist.
                    if (associated(Temp%prev)) then
                        !  Temp is in center of linklist!!!
                        if (Temp%No .EQ. OMP_bondedWallTag) then
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
                            TempH = Nodelink(OMP_bondedWallTag,tangential_force,rolling_moment,&
                            & touching,slipping,rolling,Temp,Temp%next)
                            if (associated(Temp%next)) Temp%next%prev => TempH
                            Temp%next => TempH
                            Head(I)%No = LenNode + 1
                        end if
                    else
                        !  Temp is Head of linklist!!!
                        allocate(TempH)
                        TempH = Nodelink(OMP_bondedWallTag,tangential_force,rolling_moment,&
                        & touching,slipping,rolling,Temp,Temp%next)
                        if (associated(Temp%next)) Temp%next%prev => TempH
                        Temp%next => TempH
                        Head(I)%No = LenNode + 1
                    end if
                else
                    !  memory the separation in the Hertz linklist.
                    if (associated(Temp%prev) .AND. Temp%No.EQ.OMP_bondedWallTag) then
                        !  Temp is center of linklist!!!
                        Temp%prev%next => Temp%next
                        if(associated(Temp%next)) Temp%next%prev => Temp%prev
                        Head(I)%No = LenNode - 1
                        deallocate(Temp)
                        !  When else Temp is Head of linklist!!!
                    end if
                end if
            end if
        end do
        !$OMP END PARALLEL DO
    end do
    
    do K = 1,3
        bondedWallF(K) = bondedWallF(K)/bondedWallBody
        temp_bondedWallFM(K) = bondedWallMatI(K,1)*bondedWallFM(1) + bondedWallMatI(K,2)*bondedWallFM(2) + bondedWallMatI(K,3)*bondedWallFM(3)
    end do
    bondedWallWdotB(1) = (temp_bondedWallFM(1) - (bondedWallInertia(3)-bondedWallInertia(2))*bondedWallWB(2)*bondedWallWB(3))/bondedWallInertia(1)
    bondedWallWdotB(2) = (temp_bondedWallFM(2) - (bondedWallInertia(1)-bondedWallInertia(3))*bondedWallWB(3)*bondedWallWB(1))/bondedWallInertia(2)
    bondedWallWdotB(3) = (temp_bondedWallFM(3) - (bondedWallInertia(2)-bondedWallInertia(1))*bondedWallWB(1)*bondedWallWB(2))/bondedWallInertia(3)
    
    return
    end
