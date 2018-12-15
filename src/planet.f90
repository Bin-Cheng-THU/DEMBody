    !********************************************************************
    !     DEMBody 5.0
    !     ***********
    !
    !     Planet system force.
    !     --------------------------
    !
    !     @Gravity of Pan
    !     @Gravity of Saturn
    !     @Convected inertial force
    !     @Coriolis force
    !
    !     Note that considering the movement of Pan nucleus, 
    !     the gravity of Pan must be transferred from the origin point.
    !     So we use Gravity Body instead.
    !
    !********************************************************************
    subroutine planet()

    use global
    implicit none

    integer I,J,K
    real(kind=8) rBody     !  Pan to dust
    real(kind=8) rBodyV(3) !  Pan to dust, vector
    real(kind=8) rDust     !  Saturn to dust
    real(kind=8) rDustV(3) !  Saturn to dust, vector
    real(kind=8) rP        !  Saturn to Pan Origin, scalar
    real(kind=8) center    !  Intermediate variable

    !$OMP PARALLEL DO PRIVATE(I,K,rBody,rBodyV,rDust,rDustV,rP,center)
    do I = 1,N
        !  Note that after DEMBody 4.0 we use GravBody instead.
        !!  Gravity of Pan
        !do K = 1,3
        !    rBodyV(K) = X(K,I) - X(K,N) ! whether we consider the movement of Pan
        !end do
        !rBody = sqrt(rBodyV(1)*rBodyV(1) + rBodyV(2)*rBodyV(2) + rBodyV(3)*rBodyV(3))
        !center = -muP/(rBody*rBody*rBody)
        !do K = 1,3
        !    F(K,I) = F(K,I) + center*rBodyV(K)
        !end do

        !  Gravity of Saturn
        do K = 1,3
            rDustV(K) = rOrig(K) + X(K,I)
        end do
        rDust = sqrt(rDustV(1)*rDustV(1) + rDustV(2)*rDustV(2) + rDustV(3)*rDustV(3))
        center = -muS/(rDust*rDust*rDust)
        
        do K = 1,3
            F(K,I) = F(K,I) + center*rDustV(K)
        end do

        !  Centrifugal force
        do K = 1,2
            F(K,I) = F(K,I) + omega*omega*X(K,I)
        end do

        !  Rotation transport force
        rP = sqrt(rOrig(1)*rOrig(1) + rOrig(2)*rOrig(2) + rOrig(3)*rOrig(3))
        center = muS/(rP*rP*rP)
        do K = 1,3
            F(K,I) = F(K,I) + center*rOrig(K)
        end do

        !  Coriolis force
        F(1,I) = F(1,I) + 2.0D0*omega*Xdot(2,I)
        F(2,I) = F(2,I) - 2.0D0*omega*Xdot(1,I)        
    end do
    !$OMP END PARALLEL DO        

    return
    end