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

    if (isBiDisperse) then

        !  COM of biDisperse particles
        localYORPcenter = 0.0D0
        !$OMP PARALLEL DO PRIVATE(I,K) REDUCTION(+:localYORPcenter)
        do I = 1,biDisperseNum
            do K = 1,3
                localYORPcenter(K) = localYORPcenter(K) + biDisperseX(K,I)*biDisperseBody(I)
            end do
        end do
        !$OMP END PARALLEL DO
        do K = 1,3
            localYORPcenter(K) = localYORPcenter(K)/localYORPmass
        end do

        !$OMP PARALLEL DO PRIVATE(I,K,dist)
        do I = 1,biDisperseNum
            !  Distance vector
            do K = 1,3
                dist(K) = biDisperseX(K,I) - localYORPcenter(K)
            end do

            !  Centrifugal force
            biDisperseF(1,I) = biDisperseF(1,I) + localYORPOmega*localYORPOmega*dist(1)
            biDisperseF(2,I) = biDisperseF(2,I) + localYORPOmega*localYORPOmega*dist(2)
            biDisperseF(3,I) = biDisperseF(3,I) + 0.0D0
    
            !  Coriolis force
            biDisperseF(1,I) = biDisperseF(1,I) + 2.0D0*localYORPOmega*biDisperseXdot(2,I)
            biDisperseF(2,I) = biDisperseF(2,I) - 2.0D0*localYORPOmega*biDisperseXdot(1,I)
            biDisperseF(3,I) = biDisperseF(3,I) + 0.0D0
            
            !!  Gravity
            !do K = 1,3
            !    biDisperseF(K,I) = biDisperseF(K,I) - 4.0D0/3.0D0*PI*GravConst*localYORPdensity*biDisperseX(K,I)
            !end do
        end do
        !$OMP END PARALLEL DO
    end if


    !$OMP PARALLEL DO PRIVATE(I,K,dist,distS,distL,center)
    do I = 1,N
        !  Distance vector
        do K = 1,3
            dist(K) = X(K,I) - localYORPcenter(K)
        end do
        distL = dist(1)**2 + dist(2)**2 + dist(3)**2
        distS = sqrt(distL)
        if (distS .GE. localYORPradius) then
            center = GravConst*localYORPmass/distL/distS
        else
            center = 4.0D0/3.0D0*PI*GravConst*localYORPdensity
        end if
        do K = 1,3
            F(K,I) = F(K,I) - center*dist(K)
        end do

        !  Centrifugal force
        F(1,I) = F(1,I) + localYORPOmega*localYORPOmega*dist(1)
        F(2,I) = F(2,I) + localYORPOmega*localYORPOmega*dist(2)
        F(3,I) = F(3,I) + 0.0D0
        

        !  Coriolis force
        F(1,I) = F(1,I) + 2.0D0*localYORPOmega*Xdot(2,I)
        F(2,I) = F(2,I) - 2.0D0*localYORPOmega*Xdot(1,I)
        F(3,I) = F(3,I) + 0.0D0
    end do
    !$OMP END PARALLEL DO

    return
    end