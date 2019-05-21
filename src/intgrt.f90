    !********************************************************************
    !     DEMBody 5.2
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
    
    if (isQuaternion) then
        call attitude
    end if
    
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
    
    if (isBiDisperse) then
        call intgrtBiDisperse
    end if
        
#ifdef LatticeSearch    
    if (refreshLattice) then    
        !ostart = omp_get_wtime()
        if (isPeriodic) then
            call periodic

            if (isBiDisperse) then
                call periodicBiDisperse
            end if
        end if
        !oend = omp_get_wtime()
        !write(*,*) 'periodic', (oend-ostart)
    
        !ostart = omp_get_wtime()
        call meshGenerate
        !oend = omp_get_wtime()
        !write(*,*) 'mesh', (oend-ostart)

        !ostart = omp_get_wtime()
        call latticeGenerate
        !oend = omp_get_wtime()
        !write(*,*) 'lattice', (oend-ostart)
        
        !ostart = omp_get_wtime()
        if (isBiDisperse) then
            call latticeGenerateBiDisperse
        end if    
        !oend = omp_get_wtime()
        !write(*,*) 'bidisperse lattice', (oend-ostart)
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
    
    if (isBiDisperse) then
        !$OMP PARALLEL DO PRIVATE(I,K)
        do I = 1,biDisperseNum
            do K = 1,3
                biDisperseXdot(K,I) = biDisperseXdot(K,I) + biDisperseF(K,I) * Dt /2.0D0
                biDisperseW(K,I) = biDisperseW(K,I) + biDisperseFM(K,I) * Dt /2.0D0
            end do
        end do
        !$OMP END PARALLEL DO
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