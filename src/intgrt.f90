    !********************************************************************
    !     DEMBody 2.0
    !     ***********
    !
    !     N-body integrator flow control.
    !     -------------------------------
    !     @New method to solve the intgret: mid-step velocities verlet algoritm
    !     @Calculated x & xdot. & w from the previous values
    !
    !********************************************************************
    subroutine intgrt()

    use global
    implicit none

    integer :: I,J,K
    integer :: norm

1   norm = 1

    !################         Part 1          ###################  
    !$OMP PARALLEL DO PRIVATE(I,K)
    do I = 1,N
        do K = 1,3
            X(K,I) = X(K,I) + Xdot(K,I) * Dt +  F(K,I) * Dt * Dt /2.0
            Xdot(K,I) = Xdot(K,I) + F(K,I) * Dt /2.0
            W(K,I) = W(K,I) + FM(K,I) * Dt /2.0
        end do
    end do
    !$OMP END PARALLEL DO

    if (isQuaternion) then
        call attitude
    end if
    
    if (isPeriodic) then
        call periodic
    end if

    call meshGenerate

    !################         Part 2          ################### 
    !  calculated the force from the current x & xdot.
    call force

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
    
    if (isQuaternion) then
        call attitude
    end if

    !################         Part 4          ###################
    Time = Time + Dt 

    !  check next output time.
    if (Time .GE. Tnext) then
        go to 50
    else
        go to 1
    end if


50  return

    end