    !********************************************************************
    !     DEMBody 4.3
    !     ***********
    !
    !     N-body integrator of Quaternion.
    !     -------------------------------
    !     @Refresh the quaternion of GravBody in T/2
    !
    !********************************************************************
    subroutine attitudeGravBody()

    use global
    implicit none

    integer :: J,K
    real(8) :: gravBodyQdot(4)
    real(8) :: norm
    real(8) :: gravBodyAMat_i(3,3),gravBodyWb(3)

    gravBodyAMat_i(1,1) = gravBodyQ(1)**2-gravBodyQ(2)**2-gravBodyQ(3)**2+gravBodyQ(4)**2
    gravBodyAMat_i(2,1) = 2.0D0*(gravBodyQ(1)*gravBodyQ(2)-gravBodyQ(3)*gravBodyQ(4))
    gravBodyAMat_i(3,1) = 2.0D0*(gravBodyQ(1)*gravBodyQ(3)+gravBodyQ(2)*gravBodyQ(4))
    gravBodyAMat_i(1,2) = 2.0D0*(gravBodyQ(1)*gravBodyQ(2)+gravBodyQ(3)*gravBodyQ(4))
    gravBodyAMat_i(2,2) = -gravBodyQ(1)**2+gravBodyQ(2)**2-gravBodyQ(3)**2+gravBodyQ(4)**2
    gravBodyAMat_i(3,2) = 2.0D0*(gravBodyQ(2)*gravBodyQ(3)-gravBodyQ(1)*gravBodyQ(4))
    gravBodyAMat_i(1,3) = 2.0D0*(gravBodyQ(1)*gravBodyQ(3)-gravBodyQ(2)*gravBodyQ(4))
    gravBodyAMat_i(2,3) = 2.0D0*(gravBodyQ(2)*gravBodyQ(3)+gravBodyQ(1)*gravBodyQ(4))
    gravBodyAMat_i(3,3) = -gravBodyQ(1)**2-gravBodyQ(2)**2+gravBodyQ(3)**2+gravBodyQ(4)**2
    do K = 1,3
        gravBodyWb(K) = gravBodyAMat_i(K,1)*gravBodyW(1) + gravBodyAMat_i(K,2)*gravBodyW(2) + gravBodyAMat_i(K,3)*gravBodyW(3)
    end do
    gravBodyQdot(1) = 0.5D0 * ( gravBodyWb(1)*gravBodyQ(4) - gravBodyWb(2)*gravBodyQ(3) + gravBodyWb(3)*gravBodyQ(2))
    gravBodyQdot(2) = 0.5D0 * ( gravBodyWb(1)*gravBodyQ(3) + gravBodyWb(2)*gravBodyQ(4) - gravBodyWb(3)*gravBodyQ(1))
    gravBodyQdot(3) = 0.5D0 * (-gravBodyWb(1)*gravBodyQ(2) + gravBodyWb(2)*gravBodyQ(1) + gravBodyWb(3)*gravBodyQ(4))
    gravBodyQdot(4) = 0.5D0 * (-gravBodyWb(1)*gravBodyQ(1) - gravBodyWb(2)*gravBodyQ(2) - gravBodyWb(3)*gravBodyQ(3))
    do K = 1,4
        gravBodyQ(K) = gravBodyQ(K) + gravBodyQdot(K)*0.5D0*Dt
    end do
    norm = sqrt(gravBodyQ(1)**2 + gravBodyQ(2)**2 + gravBodyQ(3)**2 + gravBodyQ(4)**2)
    do K = 1,4
        gravBodyQ(K) = gravBodyQ(K)/norm
    end do

    return
    end