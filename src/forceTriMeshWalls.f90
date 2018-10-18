    !********************************************************************
    !     DEMBody 4.3
    !     ***********
    !
    !     Force for trimesh walls.
    !     --------------------------
    !      
    !     @Using rolling friction similar to LIGGGHTS
    !     @Using Hertz-Mindlin contact model similar to Wada
    !     @No Cohesion model in walls
    !     @No Rolling model in bonded walls
    !     @Using damping model applicable for ice ball
    !
    !********************************************************************
    subroutine forceTriMeshWalls()

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
    logical :: slipping,rolling,twisting                               !  State of friction of T
    logical :: touching                                                !  State of touch of T
    integer :: I,J,K,L,II,LenNode                                      !  Iterator
    type(Nodelink),pointer :: Temp                                     !  Temporary pointer
    type(Nodelink),pointer :: TempH                                    !  Contact pointer
    integer :: OMP_trimeshWallTag                                      !  Tag for walls in OMP
    real(kind=8) OMP_trimeshWallPoint(3),OMP_trimeshWallVectorN(3)     !  Point and Normal Vector for walls in OMP
    real(kind=8) OMP_trimeshWallVectorTx(3),OMP_trimeshWallVectorTy(3) !  Tangential Vector for walls in OMP
    real(kind=8) OMP_trimeshWallLength(4)                              !  Length parameters for walls in OMP
    logical :: enterFlag
    real(kind=8)  RVb(3),RV1,RV2,RV3,RV4,RV5,RV6,RVu,RVv,RVn(3),RVt(3),norm
    real(kind=8)  edgeFlag
    
    !character(30) :: FileName
    !write(FileName,'(F10.5)') Time
    !FileName = trim(FileName)//'TriMesh.txt'
    !open(124,FILE=FileName)
    
    do II = 1,trimeshWallNum
        !  Allocate to variables in OMP (stack memory) and refresh Points & Vectors of TriMesh Walls
        OMP_trimeshWallTag = trimeshWallTag(II)
        do K = 1,3
            OMP_trimeshWallPoint(K) = trimeshWallPoint(K,II)
            OMP_trimeshWallVectorN(K) = trimeshWallVectorN(K,II)
            OMP_trimeshWallVectorTx(K) = trimeshWallVectorTx(K,II)
            OMP_trimeshWallVectorTy(K) = trimeshWallVectorTy(K,II)
        end do
        do K = 1,4
            OMP_trimeshWallLength(K) = trimeshWallLength(K,II)
        end do
        
        !  Loop over all bodies.
        !$OMP PARALLEL DO &
        !$OMP& firstprivate(OMP_trimeshWallTag,OMP_trimeshWallPoint,OMP_trimeshWallVectorN,OMP_trimeshWallVectorTx,OMP_trimeshWallVectorTy) &
        !$OMP& firstprivate(OMP_trimeshWallLength)&
        !$OMP& PRIVATE(RVb,RV1,RV2,RV3,RV4,RV5,RV6,RVu,RVv,RVn,RVt,norm,enterFlag,edgeFlag)&
        !$OMP& PRIVATE(Temp,TempH,LenNode,&
        !$OMP& I,J,K,L,II,Dist,DistS,DistL,DistR,DistU,Vrel,Vrot,Vtot,ERR,Vnor,Vtan,&
        !$OMP& normal_force,normal_forceL,tangential_force,tangential_forceL,&
        !$OMP& rolling_moment,rolling_momentL,twisting_moment,twisting_momentL,cohesive_force,Ap,An,Rij,Mij,Iij,&
        !$OMP& Kn,Cn,Ks,Cs,Kr,Cr,Kt,Ct,lnCOR,Dn,Ds,DsL,Dtheta,DthetaL,DthetaR,DthetaRL,DthetaT,DthetaTL,H,Mr,Mt,RV,&
        !$OMP& slipping,rolling,touching,twisting) SCHEDULE(DYNAMIC)
        do I = 1,N
            enterFlag = .false.
            do K = 1,3
                RV(K) = X(K,I) - OMP_trimeshWallPoint(K)
            end do
            ERR = RV(1)*OMP_trimeshWallVectorN(1) + RV(2)*OMP_trimeshWallVectorN(2) + RV(3)*OMP_trimeshWallVectorN(3)
            
            !if (I.EQ.PP) then
            !    write(124,'(F15.5,2X,2I5,2X)',advance='no') Time,I,OMP_trimeshWallTag
            !end if
            
            if (ERR.GE.-0.8D0*Dx .AND. ERR.LE.0.8D0*Dx) then
                do K = 1,3
                    RVb(K) = RV(K) - ERR*OMP_trimeshWallVectorN(K)
                end do
                RV1 = RVb(1)*OMP_trimeshWallVectorTx(1) + RVb(2)*OMP_trimeshWallVectorTx(2) + RVb(3)*OMP_trimeshWallVectorTx(3)
                RV2 = RVb(1)*OMP_trimeshWallVectorTy(1) + RVb(2)*OMP_trimeshWallVectorTy(2) + RVb(3)*OMP_trimeshWallVectorTy(3)
                RV3 = RV1 - OMP_trimeshWallLength(1)
                RV4 = RV2 - OMP_trimeshWallLength(2)
                RV5 = RV1 - OMP_trimeshWallLength(2)
                RV6 = RV2 - OMP_trimeshWallLength(3)
                RVu = (RV2*OMP_trimeshWallLength(2)-RV1*OMP_trimeshWallLength(3))/OMP_trimeshWallLength(4)
                RVv = (RV1*OMP_trimeshWallLength(2)-RV2*OMP_trimeshWallLength(1))/OMP_trimeshWallLength(4)
                if (RVu.GE.-0.8D0*Dx/sqrt(OMP_trimeshWallLength(1)) .AND. RVv.GE.-0.8D0*Dx/sqrt(OMP_trimeshWallLength(3)) &
                & .AND. ((RVu+RVv).LE.(1.0D0+0.8D0*Dx/sqrt(OMP_trimeshWallLength(1))+0.8D0*Dx/sqrt(OMP_trimeshWallLength(3)))) ) then
                    if (RVu.GE.0.0D0 .AND. RVv.GE.0.0D0 .AND. (RVu+RVv).LE.1.0D0) then
                        !  inner triangle
                        enterFlag = .true.
                        do K = 1,3
                            RVt(K) = RVb(K)
                            RVn(K) = RV(K) - RVt(K)
                        end do
                        edgeFlag = 1.0D0
                        goto 111
                    end if

                    if (RV1.LE.0.0D0 .AND. RV2.LE.0.0D0) then
                        !  adjacent angle
                        enterFlag = .true.
                        do K = 1,3
                            RVt(K) = 0.0D0
                            RVn(K) = RV(K) - RVt(K)
                        end do
                        edgeFlag = 0.2D0
                        goto 111
                    end if

                    if (RV3.GE.0.0D0 .AND. RV4.LE.RV3) then
                        !  opposite angle in Tx
                        enterFlag = .true.
                        do K = 1,3
                            RVt(K) = OMP_trimeshWallVectorTx(K)
                            RVn(K) = RV(K) - RVt(K)
                        end do
                        edgeFlag = 0.2D0
                        goto 111
                    end if

                    if (RV5.LE.RV6 .AND. RV6.GE.0.0D0) then
                        !  opposite angle in Ty
                        enterFlag = .true.
                        do K = 1,3
                            RVt(K) = OMP_trimeshWallVectorTy(K)
                            RVn(K) = RV(K) - RVt(K)
                        end do
                        edgeFlag = 0.2D0 
                        goto 111
                    end if

                    if (RV1.GT.0.0D0 .AND. RV3.LT.0.0D0 .AND. RVv.LT.0.0D0) then
                        !  adjacent edge in Tx
                        enterFlag = .true.
                        RVu = RV1/OMP_trimeshWallLength(1)
                        do K = 1,3
                            RVt(K) = RVu*OMP_trimeshWallVectorTx(K)
                            RVn(K) = RV(K) - RVt(K)
                        end do
                        edgeFlag = 0.5D0
                        goto 111
                    end if

                    if (RV2.GT.0.0D0 .AND. RV6.LT.0.0D0 .AND. RVu.LT.0.0D0) then
                        !  adjacent edge in Ty
                        enterFlag = .true.
                        RVv = RV2/OMP_trimeshWallLength(3)
                        do K = 1,3
                            RVt(K) = RVv*OMP_trimeshWallVectorTy(K)
                            RVn(K) = RV(K) - RVt(K)
                        end do
                        edgeFlag = 0.5D0
                        goto 111
                    end if

                    if (RV5.GT.RV6 .AND. RV4.GT.RV3 .AND. (RVu+RVv).GT.1.0D0) then
                        !  opposite edge
                        enterFlag = .true.
                        RVv = (RV2-RV1-OMP_trimeshWallLength(2)+OMP_trimeshWallLength(1))/(OMP_trimeshWallLength(1)+OMP_trimeshWallLength(3)-2.0*OMP_trimeshWallLength(2))
                        RVu = 1.0D0 - RVv
                        do K = 1,3
                            RVt(K) = RVu*OMP_trimeshWallVectorTx(K)+RVv*OMP_trimeshWallVectorTy(K)
                            RVn(K) = RV(K) - RVt(K)
                        end do
                        edgeFlag = 0.5D0
                        goto 111
                    end if 
                end if
            end if
            
            !if (I.EQ.PP) then
            !    write(124,'(L4,2X)',advance='no') enterFlag
            !end if
                        
111         if (enterFlag) then      
                !  normal vector
                do K = 1,3
                    Dist(K) = RVn(K)
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
                Dn = R(I) - DistL
                !  lookup this contact(or not) in the Hertz list
                LenNode = Head(I)%No
                Temp => Head(I)
                if (LenNode .NE. 0) then
                    Temp => Head(I)%next
                    do L = 1,LenNode
                        if (Temp%No .EQ. OMP_trimeshWallTag) then
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
                        else if (Temp%No.LT.OMP_trimeshWallTag .AND. associated(Temp%next)) then
                            Temp => Temp%next
                        else if (Temp%No .GT. OMP_trimeshWallTag) then
                            Temp => Temp%prev
                            exit
                        end if
                    end do
                end if
                
                !if (I.EQ.PP) then
                !    write(124,'(F15.5,2X)',advance='no') Dn
                !end if
                
                !  When collision calculate the repulsive restoring spring force which is generated along the normal and tangential according to Hertzian law
                if (Dn .GT. 0.0D0) then 

                    DistR = 1.0D0/DistL
                    !  calculate the normal vector
                    do K=1,3
                        DistU(K) = Dist(K)*DistR
                    end do
                    Ap = DistL
#ifdef HertzMindlinVisco                    
                    !  calculate material constant
                    Rij = R(I)
                    Mij = Body(I)
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
                    Rij = R(I)
                    Mij = Body(I)
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
                        Vrel(K) = Xdot(K,I)
                    end do
                    !  rotate relative velocity
                    Vrot(1) = (DistU(2)*W(3,I) - DistU(3)*W(2,I))*Ap
                    Vrot(2) = (DistU(3)*W(1,I) - DistU(1)*W(3,I))*Ap
                    Vrot(3) = (DistU(1)*W(2,I) - DistU(2)*W(1,I))*Ap
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

                    !  tangential force of Particle I
                    do K = 1,3
                        tangential_force(K) = - Ks*Ds(K) + Cs*Vtan(K) + H(K)
                    end do
                    tangential_forceL = sqrt(tangential_force(1)*tangential_force(1) + tangential_force(2)*tangential_force(2) + tangential_force(3)*tangential_force(3))

                    if (slipping) then  !  Have slipped
                        if (DsL .GT. 1.0e-8) then  !  Still slipping
                            do K = 1,3
                                tangential_force(K) = -m_mu_d*normal_forceL*Ds(K)/DsL  !  Particle I
                            end do
                        else  !  Approach sticking
                            do K = 1,3
                                tangential_force(K) = 0.0D0  !  Particle I
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
                        F(K,I) = edgeFlag*normal_force(K) + edgeFlag*tangential_force(K) + F(K,I)
                    end do
                    !  Apply moment                  
                    FM(1,I) = - edgeFlag*Ap*(DistU(2)*tangential_force(3)-DistU(3)*tangential_force(2)) + FM(1,I) 
                    FM(2,I) = - edgeFlag*Ap*(DistU(3)*tangential_force(1)-DistU(1)*tangential_force(3)) + FM(2,I) 
                    FM(3,I) = - edgeFlag*Ap*(DistU(1)*tangential_force(2)-DistU(2)*tangential_force(1)) + FM(3,I) 
                    
                    !if (I.EQ.PP) then
                    !    write(124,'(2F15.5,2X)',advance='no') normal_forceL,tangential_forceL
                    !end if
                    
                    !  rolling
                    do K = 1,3
                        Dtheta(K) = (W(K,I))*Dt
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

                    !  rolling moment of Particle I
                    do K = 1,3
                        rolling_moment(K) = - Kr*DthetaR(K) + Mr(K)
                    end do
                    rolling_momentL = sqrt(rolling_moment(1)*rolling_moment(1) + rolling_moment(2)*rolling_moment(2) + rolling_moment(3)*rolling_moment(3))                

                    if (rolling) then  !  Have rolled
                        if (DthetaRL .GT. 1.0e-8) then  !  Still slipping
                            do K = 1,3
                                rolling_moment(K) = -2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL  !  Particle I
                            end do
                        else  !  Approach sticking
                            do K = 1,3
                                rolling_moment(K) = 0.0D0  !  Particle I
                            end do
                            rolling = .false.
                        end if
                    else
                        if (rolling_momentL .GT. 2.1D0*0.25D0*m_Beta*Rij*normal_forceL) then  !  Rolling
                            rolling = .true.
                            if (DthetaRL .GT. 1.0e-14) then
                                do K = 1,3
                                    rolling_moment(K) = -2.1D0*0.25D0*m_Beta*Rij*normal_forceL*DthetaR(K)/DthetaRL
                                end do
                            else
                                do K = 1,3
                                    rolling_moment(K) = 0.0D0
                                end do
                            end if
                        else
                            rolling = .false. !  Sticking
                            do K = 1,3
                                rolling_moment(K) = rolling_moment(K) + Cr*DthetaR(K)/Dt
                            end do
                        end if
                    end if                                  

                    !  twisting moment of Particle I
                    do K = 1,3
                        twisting_moment(K) = - Kt*DthetaT(K) + Mt(K)
                    end do
                    twisting_momentL = sqrt(twisting_moment(1)*twisting_moment(1) + twisting_moment(2)*twisting_moment(2) + twisting_moment(3)*twisting_moment(3))                

                    if (twisting) then  !  Have twisted
                        if (DthetaTL .GT. 1.0e-8) then  !  Still slipping
                            do K = 1,3
                                twisting_moment(K) = -0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL  !  Particle I
                            end do
                        else  !  Approach sticking
                            do K = 1,3
                                twisting_moment(K) = 0.0D0  !  Particle I
                            end do
                            twisting = .false.
                        end if
                    else
                        if (twisting_momentL .GT. 0.65D0*m_mu_s*m_Beta*Rij*normal_forceL) then  !  Rolling
                            twisting = .true.
                            if (DthetaTL .GT. 1.0e-14) then
                                do K = 1,3
                                    twisting_moment(K) = -0.65D0*m_mu_d*m_Beta*Rij*normal_forceL*DthetaT(K)/DthetaTL
                                end do
                            else
                                do K = 1,3
                                    twisting_moment(K) = 0.0D0
                                end do
                            end if
                        else
                            twisting = .false. !  Sticking
                            do K = 1,3
                                twisting_moment(K) = twisting_moment(K) - Ct*DthetaT(K)/Dt
                            end do
                        end if
                    end if 
                    
                    !  Apply moment
                    do K = 1,3
                        FM(K,I) = edgeFlag*rolling_moment(K) + edgeFlag*twisting_moment(K) + FM(K,I)
                    end do

                    !  cohesive force
                    do K = 1,3
                        cohesive_force(K) = - m_c*(m_Beta*Rij)**2*DistU(K)
                    end do
                    do K = 1,3
                        F(K,I) = edgeFlag*cohesive_force(K) + F(K,I)
                    end do
                            
                    !  memory the contact in the Hertz linklist.
                    if (associated(Temp%prev)) then
                        !  Temp is in center of linklist!!!
                        if (Temp%No .EQ. OMP_trimeshWallTag) then
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
                            Temp%recordTime = Time + Dt
                        else
                            !  First contacted.
                            allocate(TempH)
                            TempH = Nodelink(OMP_trimeshWallTag,Time+Dt,tangential_force,rolling_moment,twisting_moment,&
                            & touching,slipping,rolling,twisting,Temp,Temp%next)
                            if (associated(Temp%next)) Temp%next%prev => TempH
                            Temp%next => TempH
                            Head(I)%No = LenNode + 1
                        end if
                    else
                        !  Temp is Head of linklist!!!
                        allocate(TempH)
                        TempH = Nodelink(OMP_trimeshWallTag,Time+Dt,tangential_force,rolling_moment,twisting_moment,&
                        & touching,slipping,rolling,twisting,Temp,Temp%next)
                        if (associated(Temp%next)) Temp%next%prev => TempH
                        Temp%next => TempH
                        Head(I)%No = LenNode + 1
                    end if
                else
                    !  memory the separation in the Hertz linklist.
                    if (associated(Temp%prev) .AND. Temp%No.EQ.OMP_trimeshWallTag) then
                        !  Temp is center of linklist!!!
                        Temp%prev%next => Temp%next
                        if(associated(Temp%next)) Temp%next%prev => Temp%prev
                        Head(I)%No = LenNode - 1
                        deallocate(Temp)
                        !  When else Temp is Head of linklist!!!
                    end if
                end if
            end if
            
            !if (I.EQ.PP) then
            !    write(124,*)
            !end if
            
        end do
        !$OMP END PARALLEL DO
    end do
    
    !write(*,*) Time,(F(K,PP),K=1,3)
    
    return
    end
