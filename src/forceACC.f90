    !********************************************************************
    !     DEMBody 6.0
    !     ***********
    !
    !     Force to Acceleration.
    !     --------------------------
    !
    !     @Add constant gravity
    !     @Force and moment normalization
    !     @Bond-force generation
    !
    !********************************************************************
    subroutine forceACC()

    use global
    implicit none

    integer I,J,K
    real(8) :: accMag
    real(8) :: t_r(3),t_moment(3)
    integer :: index

    bondF = 0.0D0
    bondFM = 0.0D0

    !$OMP PARALLEL DO PRIVATE(I,J,K,t_r)
    do I = 1,bondN
        do J = bondRange(1,I),bondRange(2,I)
            !  r
            do K = 1,3
                t_r(K) = bondMatB(K,1,I)*bondParticles(1,J) + bondMatB(K,2,I)*bondParticles(2,J) + bondMatB(K,3,I)*bondParticles(3,J)
            end do
            !  F
            do K = 1,3
                bondF(K,I) = bondF(K,I) + F(K,J)
            end do
            !  FM
            bondFM(1,I) = bondFM(1,I) + t_r(2)*F(3,J) - t_r(3)*F(2,J) + FM(1,J)
            bondFM(2,I) = bondFM(2,I) + t_r(3)*F(1,J) - t_r(1)*F(3,J) + FM(2,J)
            bondFM(3,I) = bondFM(3,I) + t_r(1)*F(2,J) - t_r(2)*F(1,J) + FM(3,J)
        end do
    end do
    !$OMP END PARALLEL DO
        
    !$OMP PARALLEL DO PRIVATE(I,K,t_moment)
    do I = 1,bondN
        do K = 1,3
            bondF(K,I) = bondF(K,I)/bondBody(I) + G(K)
            t_moment(K) = bondMatI(K,1,I)*bondFM(1,I) + bondMatI(K,2,I)*bondFM(2,I) + bondMatI(K,3,I)*bondFM(3,I)
        end do
        bondWdotB(1,I) = (t_moment(1) - (bondInertia(3,I)-bondInertia(2,I))*bondWB(2,I)*bondWB(3,I))/bondInertia(1,I)
        bondWdotB(2,I) = (t_moment(2) - (bondInertia(1,I)-bondInertia(3,I))*bondWB(3,I)*bondWB(1,I))/bondInertia(2,I)
        bondWdotB(3,I) = (t_moment(3) - (bondInertia(2,I)-bondInertia(1,I))*bondWB(1,I)*bondWB(2,I))/bondInertia(3,I)
    end do
    !$OMP END PARALLEL DO

    return
    end