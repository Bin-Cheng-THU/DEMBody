    !********************************************************************
    !     DEMBody 8.2
    !     ***********
    !
    !     YORP system force in local frame.
    !     --------------------------
    !
    !     @Convected inertial force
    !     @Coriolis force
    !     @Gravity
    !
    !     Note that here we use a simply sphere to conduct the gravity of
    !     steroids. It should be replaced with more accurate method instead.
    !
    !********************************************************************
    subroutine forceLocalYORP()

    use global
    implicit none

    integer I,J,K
    
    real(8) :: dist(3),distS,distL,center

    !$OMP PARALLEL DO PRIVATE(I,K)
    do I = 1,N
        !  Centrifugal force
        do K = 1,2
            F(K,I) = F(K,I) + localYORPOmega*localYORPOmega*X(K,I)
        end do

        !  Coriolis force
        F(1,I) = F(1,I) + 2.0D0*localYORPOmega*Xdot(2,I)
        F(2,I) = F(2,I) - 2.0D0*localYORPOmega*Xdot(1,I)
        
        !  Gravity
        do K = 1,3
            F(K,I) = F(K,I) - 4.0D0/3.0D0*PI*GravConst*localYORPdensity*X(K,I)
        end do
    end do
    !$OMP END PARALLEL DO

    if (isBiDisperse) then
        !$OMP PARALLEL DO PRIVATE(I,K)
        do I = 1,biDisperseNum
            !  Centrifugal force
            do K = 1,2
                biDisperseF(K,I) = biDisperseF(K,I) + localYORPOmega*localYORPOmega*biDisperseX(K,I)
            end do

            !  Coriolis force
            biDisperseF(1,I) = biDisperseF(1,I) + 2.0D0*localYORPOmega*biDisperseXdot(2,I)
            biDisperseF(2,I) = biDisperseF(2,I) - 2.0D0*localYORPOmega*biDisperseXdot(1,I)
            
            !  Gravity
            do K = 1,3
                biDisperseF(K,I) = biDisperseF(K,I) - 4.0D0/3.0D0*PI*GravConst*localYORPdensity*biDisperseX(K,I)
            end do
        end do
        !$OMP END PARALLEL DO
    end if

    return
    end