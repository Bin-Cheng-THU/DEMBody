    !********************************************************************
    !     DEMbody 2.0
    !     ***********
    !
    !     N-body integrator of Quaternion.
    !     -------------------------------
    !     @Refresh the quaternion of particles in T/2
    !
    !********************************************************************
    subroutine attitude()

    use global
    implicit none

    integer :: I,J,K
    real(8) :: QuaternionDot(4)
    real(8) :: norm

    !$OMP PARALLEL DO PRIVATE(I,K,norm)
    do I = 1,N
        QuaternionDot(1) = 0.5 * ( W(1,I)*Quaternion(4,I) - W(2,I)*Quaternion(3,I) + W(3,I)*Quaternion(2,I))
        QuaternionDot(2) = 0.5 * ( W(1,I)*Quaternion(3,I) + W(2,I)*Quaternion(4,I) - W(3,I)*Quaternion(1,I))
        QuaternionDot(3) = 0.5 * (-W(1,I)*Quaternion(2,I) + W(2,I)*Quaternion(1,I) + W(3,I)*Quaternion(4,I))
        QuaternionDot(4) = 0.5 * (-W(1,I)*Quaternion(1,I) - W(2,I)*Quaternion(2,I) - W(3,I)*Quaternion(3,I))
        do K = 1,4
            Quaternion(K,I) = Quaternion(K,I) + QuaternionDot(K)*0.5*Dt
        end do
        norm = sqrt(Quaternion(1,I)**2 + Quaternion(2,I)**2 + Quaternion(3,I)**2 + Quaternion(4,I)**2)
        do K = 1,4
            Quaternion(K,I) = Quaternion(K,I)/norm
        end do
    end do
    !$OMP END PARALLEL DO

    return
    end