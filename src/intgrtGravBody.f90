    !********************************************************************
    !     DEMBody 5.0
    !     ***********
    !
    !     N-body integrator of GravBody.
    !     -------------------------------
    !     @Refresh the position of GravBody in T
    !     @Refresh the velocity/angular velocity of GravBody in T/2
    !     @Refresh the quaternion of GravBody in T
    !
    !********************************************************************
    subroutine intgrtGravBody()

    use global
    implicit none

    integer :: K
    real(8) :: gravBodyQdot(4),gravBodyQdotdot(4)
    real(8) :: norm
    real(8) :: gravBodyAMat_i(3,3)
    real(8) :: gravBodyWb(3), gravBodyWdotb(3)

    !  refresh X, V
    do K = 1,3
        gravBodyX(K) = gravBodyX(K) + gravBodyXdot(K) * Dt +  gravBodyF(K) * Dt * Dt /2.0D0
        gravBodyXdot(K) = gravBodyXdot(K) + gravBodyF(K) * Dt /2.0D0
    end do
    
    !  refresh Q
    !  Q -> AMat_i
    gravBodyAMat_i(1,1) = gravBodyQ(1)**2-gravBodyQ(2)**2-gravBodyQ(3)**2+gravBodyQ(4)**2
    gravBodyAMat_i(2,1) = 2.0D0*(gravBodyQ(1)*gravBodyQ(2)-gravBodyQ(3)*gravBodyQ(4))
    gravBodyAMat_i(3,1) = 2.0D0*(gravBodyQ(1)*gravBodyQ(3)+gravBodyQ(2)*gravBodyQ(4))
    gravBodyAMat_i(1,2) = 2.0D0*(gravBodyQ(1)*gravBodyQ(2)+gravBodyQ(3)*gravBodyQ(4))
    gravBodyAMat_i(2,2) = -gravBodyQ(1)**2+gravBodyQ(2)**2-gravBodyQ(3)**2+gravBodyQ(4)**2
    gravBodyAMat_i(3,2) = 2.0D0*(gravBodyQ(2)*gravBodyQ(3)-gravBodyQ(1)*gravBodyQ(4))
    gravBodyAMat_i(1,3) = 2.0D0*(gravBodyQ(1)*gravBodyQ(3)-gravBodyQ(2)*gravBodyQ(4))
    gravBodyAMat_i(2,3) = 2.0D0*(gravBodyQ(2)*gravBodyQ(3)+gravBodyQ(1)*gravBodyQ(4))
    gravBodyAMat_i(3,3) = -gravBodyQ(1)**2-gravBodyQ(2)**2+gravBodyQ(3)**2+gravBodyQ(4)**2
    !  omega_i -> omega_b
    do K = 1,3
        gravBodyWb(K) = gravBodyAMat_i(K,1)*gravBodyW(1) + gravBodyAMat_i(K,2)*gravBodyW(2) + gravBodyAMat_i(K,3)*gravBodyW(3)
    end do
    !  omegadot_i -> omegadot_b
    do K = 1,3
        gravBodyWdotb(K) = gravBodyAMat_i(K,1)*gravBodyFM(1) + gravBodyAMat_i(K,2)*gravBodyFM(2) + gravBodyAMat_i(K,3)*gravBodyFM(3)
    end do
    !  Q, omega -> Qdot
    gravBodyQdot(1) = 0.5D0 * ( gravBodyWb(1)*gravBodyQ(4) - gravBodyWb(2)*gravBodyQ(3) + gravBodyWb(3)*gravBodyQ(2))
    gravBodyQdot(2) = 0.5D0 * ( gravBodyWb(1)*gravBodyQ(3) + gravBodyWb(2)*gravBodyQ(4) - gravBodyWb(3)*gravBodyQ(1))
    gravBodyQdot(3) = 0.5D0 * (-gravBodyWb(1)*gravBodyQ(2) + gravBodyWb(2)*gravBodyQ(1) + gravBodyWb(3)*gravBodyQ(4))
    gravBodyQdot(4) = 0.5D0 * (-gravBodyWb(1)*gravBodyQ(1) - gravBodyWb(2)*gravBodyQ(2) - gravBodyWb(3)*gravBodyQ(3))
    !  Q, omegadot, Qdot -> Qdotdot
    gravBodyQdotdot(1) = ( gravBodyWdotb(1)*gravBodyQ(4) - gravBodyWdotb(2)*gravBodyQ(3) + gravBodyWdotb(3)*gravBodyQ(2))*0.5D0 + ( gravBodyWb(1)*gravBodyQdot(4) - gravBodyWb(2)*gravBodyQdot(3) + gravBodyWb(3)*gravBodyQdot(2))*0.5D0
    gravBodyQdotdot(2) = ( gravBodyWdotb(1)*gravBodyQ(3) + gravBodyWdotb(2)*gravBodyQ(4) - gravBodyWdotb(3)*gravBodyQ(1))*0.5D0 + ( gravBodyWb(1)*gravBodyQdot(3) + gravBodyWb(2)*gravBodyQdot(4) - gravBodyWb(3)*gravBodyQdot(1))*0.5D0
    gravBodyQdotdot(3) = (-gravBodyWdotb(1)*gravBodyQ(2) + gravBodyWdotb(2)*gravBodyQ(1) + gravBodyWdotb(3)*gravBodyQ(4))*0.5D0 + (-gravBodyWb(1)*gravBodyQdot(2) + gravBodyWb(2)*gravBodyQdot(1) + gravBodyWb(3)*gravBodyQdot(4))*0.5D0
    gravBodyQdotdot(4) = (-gravBodyWdotb(1)*gravBodyQ(1) - gravBodyWdotb(2)*gravBodyQ(2) - gravBodyWdotb(3)*gravBodyQ(3))*0.5D0 + (-gravBodyWb(1)*gravBodyQdot(1) - gravBodyWb(2)*gravBodyQdot(2) - gravBodyWb(3)*gravBodyQdot(3))*0.5D0
    do K = 1,4
        gravBodyQ(K) = gravBodyQ(K) + gravBodyQdot(K) * Dt + gravBodyQdotdot(K) * Dt * Dt /2.0D0
    end do
    norm = sqrt(gravBodyQ(1)**2 + gravBodyQ(2)**2 + gravBodyQ(3)**2 + gravBodyQ(4)**2)
    do K = 1,4
        gravBodyQ(K) = gravBodyQ(K)/norm
    end do
    
    !  refresh omega
    do K = 1,3
        gravBodyW(K) = gravBodyW(K) + gravBodyFM(K) * Dt /2.0D0
    end do

    return
    end