    !********************************************************************
    !     DEMBody 4.0
    !     ***********
    !
    !     N-body integrator flow control.
    !     -------------------------------
    !     @New method to solve the intgret: mid-step velocities verlet algoritm
    !     @Calculated x & xdot & w from the previous values.

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

1   norm = 1
    
    !o1 = omp_get_wtime()

    !################         Part 1          ###################  
    !  calculated the force from the current x & xdot.
    !$OMP PARALLEL DO PRIVATE(I,K)
    do I = 1,N
        do K = 1,3
            X(K,I) = X(K,I) + Xdot(K,I) * Dt +  F(K,I) * Dt * Dt /2.0
            Xdot(K,I) = Xdot(K,I) + F(K,I) * Dt /2.0
            W(K,I) = W(K,I) + FM(K,I) * Dt /2.0
        end do
    end do
    !$OMP END PARALLEL DO
    
    if (isBondedWall) then
        call attitudeBondedWalls
    end if
    
    if (isGravBody) then
        do K = 1,3
            gravBodyX(K) = gravBodyX(K) + gravBodyXdot(K) * Dt +  gravBodyF(K) * Dt * Dt /2.0
            gravBodyXdot(K) = gravBodyXdot(K) + gravBodyF(K) * Dt /2.0
            gravBodyW(K) = gravBodyW(K) + gravBodyFM(K) * Dt /2.0
        end do
        call attitudeGravBody
    end if

    if (isQuaternion) then
        call attitude
    end if
    
    if (isPeriodic) then
        call periodic
    end if

    !ostart = omp_get_wtime()
    call meshGenerate
    !oend = omp_get_wtime()
    !write(*,*) 'mesh', (oend-ostart)
    
    !ostart = omp_get_wtime()
    call latticeGenerate
    !oend = omp_get_wtime()
    !write(*,*) 'lattice', (oend-ostart)

    !################         Part 2          ################### 
    !  calculated the force from the current x & xdot.
    call forceLattice

    !################         Part 3          ################### 
    !  calculated the current xdot
    !$OMP PARALLEL DO PRIVATE(I,K)
    do I = 1,N
        do K = 1,3
            Xdot(K,I) = Xdot(K,I) + F(K,I) * Dt /2.0
            W(K,I) = W(K,I) + FM(K,I) * Dt /2.0            
        end do       
    end do
    !$OMP END PARALLEL DO
    
    if (isBondedWall) then
        do K = 1,3
            bondedWallXdot(K) = bondedWallXdot(K) + bondedWallF(K) * Dt /2.0
            bondedWallWB(K) = bondedWallWB(K) + bondedWallWdotB(K) * Dt /2.0
        end do
        do K = 1,3
            bondedWallW(K) = bondedWallMatB(K,1)*bondedWallWB(1) + bondedWallMatB(K,2)*bondedWallWB(2) + bondedWallMatB(K,3)*bondedWallWB(3)
        end do
    end if
    
    if (isGravBody) then
        do K = 1,3
            gravBodyXdot(K) = gravBodyXdot(K) + gravBodyF(K) * Dt /2.0
            gravBodyW(K) = gravBodyW(K) + gravBodyFM(K) * Dt /2.0
        end do
        call attitudeGravBody
    end if
    
    if (isQuaternion) then
        call attitude
    end if

    !################         Part 4          ###################
    Time = Time + Dt 
    
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