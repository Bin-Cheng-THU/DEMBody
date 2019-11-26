    !********************************************************************
    !     DEMBody 7.0
    !     ***********
    !
    !     N-body integrator flow control.
    !     -------------------------------
    !     @New method to solve the intgret: mid-step velocities verlet algoritm
    !     @Calculated x & xdot & w & Q from the previous values.
    !
    !********************************************************************
    subroutine intgrt()

    use global
    use omp_lib
    implicit none

    integer :: I,J,K
    integer :: norm
    real(8) :: ostart,oend
    real(8) :: o1,o2
    real(8) :: dist(3),distL
    real(8) :: Tmoving,Dmoving,Vmoving

1   norm = 1
    
    !o1 = omp_get_wtime()
    
    !################         Part 1          ###################
    Time = Time + Dt 
    currentStep = currentStep + 1

    !################         Part 2          ###################  
    !  calculated the force from the current x & xdot.
    !$OMP PARALLEL DO PRIVATE(I,K,dist,distL)
    do I = 1,N
        do K = 1,3
            dist(K) = Xdot(K,I) * Dt +  F(K,I) * Dt * Dt /2.0D0
            XT(K,I) = XT(K,I) + dist(K)            
            X(K,I) = X(K,I) + dist(K)
            Xdot(K,I) = Xdot(K,I) + F(K,I) * Dt /2.0D0
        end do
#ifdef LatticeSearch        
        distL = XT(1,I)*XT(1,I) + XT(2,I)*XT(2,I) + XT(3,I)*XT(3,I)
        if (distL .GE. verlet) then
            refreshLattice = .true.
        end if
#endif LatticeSearch        
    end do
    !$OMP END PARALLEL DO
    
    !  refresh quaternion if needed
    if (isQuaternion) then
        call attitude
    end if

    !$OMP PARALLEL DO PRIVATE(I,K)
    do I = 1,N
        !  refresh omega
        do K = 1,3
            W(K,I) = W(K,I) + FM(K,I) * Dt /2.0D0
        end do
    end do
    !$OMP END PARALLEL DO
    
    if (isMovingWall) then
        Tmoving = Time - movingWallTstart
        if (Tmoving.GE.0.0D0 .AND. Tmoving.LT.(movingWallTend-movingWallTstart)) then
            Dmoving = movingWallA*sin(movingWallOmega*Tmoving)
            Vmoving = movingWallA*movingWallOmega*cos(movingWallOmega*Tmoving)
        else if (Tmoving.GE.(movingWallTend-movingWallTstart)) then
            Dmoving = movingWallA*sin(movingWallOmega*(movingWallTend-movingWallTstart))
            Vmoving = 0.0D0
        else
            Dmoving = 0.0D0
            Vmoving = 0.0D0
        end if
        do I = 1,movingWallNum
            do K = 1,3
                movingWallPoint(K,I) = movingWallPointInit(K,I) + Dmoving*movingWallNormal(K)
                movingWallVelocity(K,I) = Vmoving*movingWallNormal(K)
            end do
        end do
        !do I = 1,movingWallNum
        !    do K = 1,3
        !        movingWallPoint(K,I) = movingWallPointStore(K,I,(currentStep+1))
        !        movingWallVector(K,I) = movingWallVectorStore(K,I,(currentStep+1))
        !        movingWallVelocity(K,I) = movingWallVelocityStore(K,I,(currentStep+1))
        !    end do
        !end do
    end if
    
    if (isBondedWall) then
        call intgrtBondedWalls
    end if
    
    if (isGravBody) then
        call intgrtGravBody
    end if
    
    if (isSphereBody) then
        call intgrtSphereBody
    end if

    if (isStretch) then
        Tmoving = Time - stretchTstart
        stretchVelX = 0.0D0
        stretchVelY = 0.0D0
        if (Tmoving.GE.0.0D0 .AND. Tmoving.LT.(stretchTend-stretchTstart)) then
            !  refresh peirodic boundary
            stretchVelX = stretchVelXInit
            stretchVelY = stretchVelYInit
            PlaSx1p(1) = PlaSx1p(1) + stretchVelX*Dt
            PlaSx2p(1) = PlaSx2p(1) - stretchVelX*Dt
            PlaSy1p(2) = PlaSy1p(2) + stretchVelY*Dt
            PlaSy2p(2) = PlaSy2p(2) - stretchVelY*Dt
            LenBoxX = abs(PlaSx1p(1) - PlaSx2p(1))
            LenBoxY = abs(PlaSy1p(2) - PlaSy2p(2))
            if (isPeriodic) then 
                !  refresh Lattice mesh
                LatMx = LatMx + stretchVelX*Dt
                LatMy = LatMy + stretchVelY*Dt
                LatDx = LenBoxX/LatNx
                LatDy = LenBoxY/LatNy
            end if
        end if
    end if

#ifdef LatticeSearch    
    if (refreshLattice) then    
        !ostart = omp_get_wtime()
        if (isPeriodic) then
            call periodic
        end if
        !oend = omp_get_wtime()
        !write(*,*) 'periodic', (oend-ostart)

        !ostart = omp_get_wtime()
        if (isMirror) then
            call mirror
        end if
        !oend = omp_get_wtime()
        !write(*,*) 'mirror', (oend-ostart)
    
        !ostart = omp_get_wtime()
        call meshGenerate
        !oend = omp_get_wtime()
        !write(*,*) 'mesh', (oend-ostart)

        !ostart = omp_get_wtime()
        call latticeGenerate
        !oend = omp_get_wtime()
        !write(*,*) 'lattice', (oend-ostart)
    end if
#endif
    
    !################         Part 3          ################### 
    !  calculated the force from the current x & xdot.
    !ostart = omp_get_wtime()
    call force
    !oend = omp_get_wtime()
    !write(*,*) 'force', (oend-ostart)

    !################         Part 4          ################### 
    !  calculated the current xdot
    !$OMP PARALLEL DO PRIVATE(I,K)
    do I = 1,N
        do K = 1,3
            Xdot(K,I) = Xdot(K,I) + F(K,I) * Dt /2.0D0
            W(K,I) = W(K,I) + FM(K,I) * Dt /2.0D0     
        end do       
    end do
    !$OMP END PARALLEL DO
    
    if (isBondedWall) then
        do K = 1,3
            bondedWallXdot(K) = bondedWallXdot(K) + bondedWallF(K) * Dt /2.0D0
            bondedWallWB(K) = bondedWallWB(K) + bondedWallWdotB(K) * Dt /2.0D0
        end do
        do K = 1,3
            bondedWallW(K) = bondedWallMatB(K,1)*bondedWallWB(1) + bondedWallMatB(K,2)*bondedWallWB(2) + bondedWallMatB(K,3)*bondedWallWB(3)
        end do
    end if
    
    if (isGravBody) then
        do K = 1,3
            gravBodyXdot(K) = gravBodyXdot(K) + gravBodyF(K) * Dt /2.0D0
            gravBodyW(K) = gravBodyW(K) + gravBodyFM(K) * Dt /2.0D0
        end do
    end if
    
    if (isSphereBody) then
        do I = 1,sphereBodyNum
            do K = 1,3
                sphereBodyXdot(K,I) = sphereBodyXdot(K,I) + sphereBodyF(K,I) * Dt /2.0D0
                sphereBodyW(K,I) = sphereBodyW(K,I) + sphereBodyFM(K,I) * Dt /2.0D0
            end do
        end do
    end if
    
    !o2 = omp_get_wtime()
    !write(*,*) 'intgrt',(o2-o1)
    !write(*,*)
    !pause

    !  check next output time.
    if (Time .GE. Tnext) then
        go to 50
    else
        go to 1
    end if


50  return

    end