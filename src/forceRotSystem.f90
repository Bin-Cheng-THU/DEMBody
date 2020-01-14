    !********************************************************************
    !     DEMBody 5.0
    !     ***********
    !
    !     Rotary system force.
    !     --------------------------
    !
    !     @Convected inertial force
    !     @Coriolis force
    !
    !
    !********************************************************************
    subroutine forceRotSystem()

    use global
    implicit none

    integer I,J,K
    
    real(8) :: dist(3),distS,distL,center

    !$OMP PARALLEL DO PRIVATE(I,K,dist,distS,distL,center)
    do I = 1,N
        !  Centrifugal force
        do K = 1,2
            F(K,I) = F(K,I) + sysOmega*sysOmega*X(K,I)
        end do

        !  Coriolis force
        F(1,I) = F(1,I) + 2.0D0*sysOmega*Xdot(2,I)
        F(2,I) = F(2,I) - 2.0D0*sysOmega*Xdot(1,I)
        
        !!  Gravity
        !do K = 1,3
        !    dist(K) = 0.0D0 - X(K,I)
        !end do
        !distS = Dist(1)*Dist(1) + Dist(2)*Dist(2) + Dist(3)*Dist(3)
        !distL = sqrt(DistS)
        !center = sysGrav/distS/distL
        !do K = 1,3
        !    F(K,I) = F(K,I) + center*dist(K)
        !end do
    end do
    !$OMP END PARALLEL DO

    if (isBiDisperse) then
        !$OMP PARALLEL DO PRIVATE(I,K,distS,distL,center)
        do I = 1,biDisperseNum
            !  Centrifugal force
            do K = 1,2
                biDisperseF(K,I) = biDisperseF(K,I) + sysOmega*sysOmega*biDisperseX(K,I)
            end do

            !  Coriolis force
            biDisperseF(1,I) = biDisperseF(1,I) + 2.0D0*sysOmega*biDisperseXdot(2,I)
            biDisperseF(2,I) = biDisperseF(2,I) - 2.0D0*sysOmega*biDisperseXdot(1,I)
            
            !!  Gravity
            !do K = 1,3
            !    dist(K) = 0.0D0 - biDisperseX(K,I)
            !end do
            !distS = Dist(1)*Dist(1) + Dist(2)*Dist(2) + Dist(3)*Dist(3)
            !distL = sqrt(DistS)
            !center = sysGrav/distS/distL
            !do K = 1,3
            !    biDisperseF(K,I) = biDisperseF(K,I) + center*dist(K)
            !end do
        end do
        !$OMP END PARALLEL DO
    end if

    return
    end