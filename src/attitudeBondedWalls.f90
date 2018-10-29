    !********************************************************************
    !     DEMBody 4.5
    !     ***********
    !
    !     N-body integrator of Quaternion.
    !     -------------------------------
    !     @Refresh the quaternion of bonded walls in T/2
    !
    !********************************************************************
    subroutine attitudeBondedWalls()
    
    use global
    implicit none
    
    real(8) :: Qdot(4),Qdotdot(4)
    real(8) :: norm
    integer :: K,II
    
    !  refresh X, V
    do K = 1,3
        bondedWallX(K) = bondedWallX(K) + bondedWallXdot(K) * Dt + bondedWallF(K) * Dt * Dt /2.0D0
        bondedWallXdot(K) = bondedWallXdot(K) + bondedWallF(K) * Dt /2.0D0
    end do
    
    !  refresh Quaternion dot
    Qdot(1) = ( bondedWallWB(1)*bondedWallQ(4) - bondedWallWB(2)*bondedWallQ(3) + bondedWallWB(3)*bondedWallQ(2))*0.5D0
    Qdot(2) = ( bondedWallWB(1)*bondedWallQ(3) + bondedWallWB(2)*bondedWallQ(4) - bondedWallWB(3)*bondedWallQ(1))*0.5D0
    Qdot(3) = (-bondedWallWB(1)*bondedWallQ(2) + bondedWallWB(2)*bondedWallQ(1) + bondedWallWB(3)*bondedWallQ(4))*0.5D0
    Qdot(4) = (-bondedWallWB(1)*bondedWallQ(1) - bondedWallWB(2)*bondedWallQ(2) - bondedWallWB(3)*bondedWallQ(3))*0.5D0
    !  refresh Quaternion dot dot
    Qdotdot(1) = ( bondedWallWdotB(1)*bondedWallQ(4) - bondedWallWdotB(2)*bondedWallQ(3) + bondedWallWdotB(3)*bondedWallQ(2))*0.5D0 + ( bondedWallWB(1)*Qdot(4) - bondedWallWB(2)*Qdot(3) + bondedWallWB(3)*Qdot(2))*0.5D0
    Qdotdot(2) = ( bondedWallWdotB(1)*bondedWallQ(3) + bondedWallWdotB(2)*bondedWallQ(4) - bondedWallWdotB(3)*bondedWallQ(1))*0.5D0 + ( bondedWallWB(1)*Qdot(3) + bondedWallWB(2)*Qdot(4) - bondedWallWB(3)*Qdot(1))*0.5D0
    Qdotdot(3) = (-bondedWallWdotB(1)*bondedWallQ(2) + bondedWallWdotB(2)*bondedWallQ(1) + bondedWallWdotB(3)*bondedWallQ(4))*0.5D0 + (-bondedWallWB(1)*Qdot(2) + bondedWallWB(2)*Qdot(1) + bondedWallWB(3)*Qdot(4))*0.5D0
    Qdotdot(4) = (-bondedWallWdotB(1)*bondedWallQ(1) - bondedWallWdotB(2)*bondedWallQ(2) - bondedWallWdotB(3)*bondedWallQ(3))*0.5D0 + (-bondedWallWB(1)*Qdot(1) - bondedWallWB(2)*Qdot(2) - bondedWallWB(3)*Qdot(3))*0.5D0
    !  refresh Quaternion
    do K = 1,4
        bondedWallQ(K) = bondedWallQ(K) + Qdot(K) * Dt + Qdotdot(K) * Dt * Dt /2.0D0
    end do
    !  Normalize Quaternion
    norm = sqrt(bondedWallQ(1)**2 + bondedWallQ(2)**2 + bondedWallQ(3)**2 + bondedWallQ(4)**2)
    do K = 1,4
        bondedWallQ(K) = bondedWallQ(K)/norm
    end do
    
    call attitudeBondedWallsQM
    
    !  refresh Angular velocity
    do K = 1,3
        bondedWallWB(K) = bondedWallWB(K) + bondedWallWdotB(K) * Dt /2.0D0
    end do
    !  refresh inertial bondwall angular velocity
    do K = 1,3
        bondedWallW(K) = bondedWallMatB(K,1)*bondedWallWB(1) + bondedWallMatB(K,2)*bondedWallWB(2) + bondedWallMatB(K,3)*bondedWallWB(3)
    end do
            
    end subroutine