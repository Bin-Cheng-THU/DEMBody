    !********************************************************************
    !     DEMBody 5.0
    !     ***********
    !
    !     N-body integrator of SphereBody.
    !     -------------------------------
    !     @Refresh the position of SphereBody in T
    !     @Refresh the velocity/angular velocity of SphereBody in T/2
    !     @Refresh the quaternion of SphereBody in T
    !
    !********************************************************************
    subroutine intgrtSphereBody()

    use global
    implicit none

    integer :: I,K
    real(8) :: sphereBodyQdot(4),sphereBodyQdotdot(4)
    real(8) :: norm
    real(8) :: sphereBodyAMat_i(3,3)
    real(8) :: sphereBodyWb(3), sphereBodyWdotb(3)

    !  refresh X, V
    do I = 1,sphereBodyNum
        do K = 1,3
            sphereBodyX(K,I) = sphereBodyX(K,I) + sphereBodyXdot(K,I) * Dt +  sphereBodyF(K,I) * Dt * Dt /2.0D0
            sphereBodyXdot(K,I) = sphereBodyXdot(K,I) + sphereBodyF(K,I) * Dt /2.0D0
        end do
    end do
    
    do I = 1,sphereBodyNum
        !  refresh Q
        !  Q -> AMat_i
        sphereBodyAMat_i(1,1) = sphereBodyQ(1,I)**2-sphereBodyQ(2,I)**2-sphereBodyQ(3,I)**2+sphereBodyQ(4,I)**2
        sphereBodyAMat_i(2,1) = 2.0D0*(sphereBodyQ(1,I)*sphereBodyQ(2,I)-sphereBodyQ(3,I)*sphereBodyQ(4,I))
        sphereBodyAMat_i(3,1) = 2.0D0*(sphereBodyQ(1,I)*sphereBodyQ(3,I)+sphereBodyQ(2,I)*sphereBodyQ(4,I))
        sphereBodyAMat_i(1,2) = 2.0D0*(sphereBodyQ(1,I)*sphereBodyQ(2,I)+sphereBodyQ(3,I)*sphereBodyQ(4,I))
        sphereBodyAMat_i(2,2) = -sphereBodyQ(1,I)**2+sphereBodyQ(2,I)**2-sphereBodyQ(3,I)**2+sphereBodyQ(4,I)**2
        sphereBodyAMat_i(3,2) = 2.0D0*(sphereBodyQ(2,I)*sphereBodyQ(3,I)-sphereBodyQ(1,I)*sphereBodyQ(4,I))
        sphereBodyAMat_i(1,3) = 2.0D0*(sphereBodyQ(1,I)*sphereBodyQ(3,I)-sphereBodyQ(2,I)*sphereBodyQ(4,I))
        sphereBodyAMat_i(2,3) = 2.0D0*(sphereBodyQ(2,I)*sphereBodyQ(3,I)+sphereBodyQ(1,I)*sphereBodyQ(4,I))
        sphereBodyAMat_i(3,3) = -sphereBodyQ(1,I)**2-sphereBodyQ(2,I)**2+sphereBodyQ(3,I)**2+sphereBodyQ(4,I)**2
        !  omega_i -> omega_b
        do K = 1,3
            sphereBodyWb(K) = sphereBodyAMat_i(K,1)*sphereBodyW(1,I) + sphereBodyAMat_i(K,2)*sphereBodyW(2,I) + sphereBodyAMat_i(K,3)*sphereBodyW(3,I)
        end do
        !  omegadot_i -> omegadot_b
        do K = 1,3
            sphereBodyWdotb(K) = sphereBodyAMat_i(K,1)*sphereBodyFM(1,I) + sphereBodyAMat_i(K,2)*sphereBodyFM(2,I) + sphereBodyAMat_i(K,3)*sphereBodyFM(3,I)
        end do
        !  Q, omega -> Qdot
        sphereBodyQdot(1) = 0.5D0 * ( sphereBodyWb(1)*sphereBodyQ(4,I) - sphereBodyWb(2)*sphereBodyQ(3,I) + sphereBodyWb(3)*sphereBodyQ(2,I))
        sphereBodyQdot(2) = 0.5D0 * ( sphereBodyWb(1)*sphereBodyQ(3,I) + sphereBodyWb(2)*sphereBodyQ(4,I) - sphereBodyWb(3)*sphereBodyQ(1,I))
        sphereBodyQdot(3) = 0.5D0 * (-sphereBodyWb(1)*sphereBodyQ(2,I) + sphereBodyWb(2)*sphereBodyQ(1,I) + sphereBodyWb(3)*sphereBodyQ(4,I))
        sphereBodyQdot(4) = 0.5D0 * (-sphereBodyWb(1)*sphereBodyQ(1,I) - sphereBodyWb(2)*sphereBodyQ(2,I) - sphereBodyWb(3)*sphereBodyQ(3,I))
        !  Q, omegadot, Qdot -> Qdotdot
        sphereBodyQdotdot(1) = ( sphereBodyWdotb(1)*sphereBodyQ(4,I) - sphereBodyWdotb(2)*sphereBodyQ(3,I) + sphereBodyWdotb(3)*sphereBodyQ(2,I))*0.5D0 + ( sphereBodyWb(1)*sphereBodyQdot(4) - sphereBodyWb(2)*sphereBodyQdot(3) + sphereBodyWb(3)*sphereBodyQdot(2))*0.5D0
        sphereBodyQdotdot(2) = ( sphereBodyWdotb(1)*sphereBodyQ(3,I) + sphereBodyWdotb(2)*sphereBodyQ(4,I) - sphereBodyWdotb(3)*sphereBodyQ(1,I))*0.5D0 + ( sphereBodyWb(1)*sphereBodyQdot(3) + sphereBodyWb(2)*sphereBodyQdot(4) - sphereBodyWb(3)*sphereBodyQdot(1))*0.5D0
        sphereBodyQdotdot(3) = (-sphereBodyWdotb(1)*sphereBodyQ(2,I) + sphereBodyWdotb(2)*sphereBodyQ(1,I) + sphereBodyWdotb(3)*sphereBodyQ(4,I))*0.5D0 + (-sphereBodyWb(1)*sphereBodyQdot(2) + sphereBodyWb(2)*sphereBodyQdot(1) + sphereBodyWb(3)*sphereBodyQdot(4))*0.5D0
        sphereBodyQdotdot(4) = (-sphereBodyWdotb(1)*sphereBodyQ(1,I) - sphereBodyWdotb(2)*sphereBodyQ(2,I) - sphereBodyWdotb(3)*sphereBodyQ(3,I))*0.5D0 + (-sphereBodyWb(1)*sphereBodyQdot(1) - sphereBodyWb(2)*sphereBodyQdot(2) - sphereBodyWb(3)*sphereBodyQdot(3))*0.5D0
        do K = 1,4
            sphereBodyQ(K,I) = sphereBodyQ(K,I) + sphereBodyQdot(K) * Dt + sphereBodyQdotdot(K) * Dt * Dt /2.0D0
        end do
        norm = sqrt(sphereBodyQ(1,I)**2 + sphereBodyQ(2,I)**2 + sphereBodyQ(3,I)**2 + sphereBodyQ(4,I)**2)
        do K = 1,4
            sphereBodyQ(K,I) = sphereBodyQ(K,I)/norm
        end do

        !  refresh omega
        do K = 1,3
            sphereBodyW(K,I) = sphereBodyW(K,I) + sphereBodyFM(K,I) * Dt /2.0D0
        end do

    end do

    return
    end