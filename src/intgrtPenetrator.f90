    !********************************************************************
    !     DEMBody 6.4
    !     ***********
    !
    !     N-body integrator of penetrator.
    !     -------------------------------
    !     @Refresh the position of bonded walls in T
    !     @Refresh the quaternion of bonded walls in T
    !     @Refresh the velocity/anglar velocity of bonded walls in T/2
    !
    !********************************************************************
    subroutine intgrtPenetrator()
    
    use global
    implicit none
    
    real(8) :: Qdot(4),Qdotdot(4)
    real(8) :: norm
    integer :: K
    
    !  refresh X, V
    do K = 1,3
        peneX(K) = peneX(K) + peneXdot(K) * Dt + peneF(K) * Dt * Dt /2.0D0
        peneXdot(K) = peneXdot(K) + peneF(K) * Dt /2.0D0
    end do
    
    !  refresh Q
    !  Q, omega -> Qdot
    Qdot(1) = ( peneWB(1)*peneQ(4) - peneWB(2)*peneQ(3) + peneWB(3)*peneQ(2))*0.5D0
    Qdot(2) = ( peneWB(1)*peneQ(3) + peneWB(2)*peneQ(4) - peneWB(3)*peneQ(1))*0.5D0
    Qdot(3) = (-peneWB(1)*peneQ(2) + peneWB(2)*peneQ(1) + peneWB(3)*peneQ(4))*0.5D0
    Qdot(4) = (-peneWB(1)*peneQ(1) - peneWB(2)*peneQ(2) - peneWB(3)*peneQ(3))*0.5D0
    !  Q, omegaDot, Qdot -> Qdotdot
    Qdotdot(1) = ( peneWdotB(1)*peneQ(4) - peneWdotB(2)*peneQ(3) + peneWdotB(3)*peneQ(2))*0.5D0 + ( peneWB(1)*Qdot(4) - peneWB(2)*Qdot(3) + peneWB(3)*Qdot(2))*0.5D0
    Qdotdot(2) = ( peneWdotB(1)*peneQ(3) + peneWdotB(2)*peneQ(4) - peneWdotB(3)*peneQ(1))*0.5D0 + ( peneWB(1)*Qdot(3) + peneWB(2)*Qdot(4) - peneWB(3)*Qdot(1))*0.5D0
    Qdotdot(3) = (-peneWdotB(1)*peneQ(2) + peneWdotB(2)*peneQ(1) + peneWdotB(3)*peneQ(4))*0.5D0 + (-peneWB(1)*Qdot(2) + peneWB(2)*Qdot(1) + peneWB(3)*Qdot(4))*0.5D0
    Qdotdot(4) = (-peneWdotB(1)*peneQ(1) - peneWdotB(2)*peneQ(2) - peneWdotB(3)*peneQ(3))*0.5D0 + (-peneWB(1)*Qdot(1) - peneWB(2)*Qdot(2) - peneWB(3)*Qdot(3))*0.5D0
    !  refresh Quaternion
    do K = 1,4
        peneQ(K) = peneQ(K) + Qdot(K) * Dt + Qdotdot(K) * Dt * Dt /2.0D0
    end do
    !  Normalize Quaternion
    norm = sqrt(peneQ(1)**2 + peneQ(2)**2 + peneQ(3)**2 + peneQ(4)**2)
    do K = 1,4
        peneQ(K) = peneQ(K)/norm
    end do
    
    call attitudeQ2M(peneQ,peneMatI,peneMatB)
    
    !  refresh Angular velocity
    do K = 1,3
        peneWB(K) = peneWB(K) + peneWdotB(K) * Dt /2.0D0
        peneVector(K) = peneMatB(K,3)
    end do
    !  refresh Inertial bondwall angular velocity
    do K = 1,3
        peneW(K) = peneMatB(K,1)*peneWB(1) + peneMatB(K,2)*peneWB(2) + peneMatB(K,3)*peneWB(3)
    end do
            
    end subroutine