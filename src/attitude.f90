    !********************************************************************
    !     DEMBody 5.0
    !     ***********
    !
    !     N-body integrator of Quaternion.
    !     -------------------------------
    !     @Refresh the quaternion of particles in T
    !     @Refresh the angular velocity of particles in T/2
    !
    !********************************************************************
    subroutine attitude()

    use global
    implicit none

    integer :: I,K
    real(8) :: QuaternionDot(4),QuaternionDotDot(4)
    real(8) :: norm
    real(8) :: OMP_AMat_i(3,3),OMP_Wb(3),OMP_Wdotb(3)

    !$OMP PARALLEL DO PRIVATE(I,K,QuaternionDot,QuaternionDotDot,norm,OMP_AMat_i,OMP_Wb,OMP_Wdotb)
    do I = 1,N
        !  refresh Q
        !  Q -> AMat_i
        OMP_AMat_i(1,1) = Quaternion(1,I)**2-Quaternion(2,I)**2-Quaternion(3,I)**2+Quaternion(4,I)**2
        OMP_AMat_i(2,1) = 2.0D0*(Quaternion(1,I)*Quaternion(2,I)-Quaternion(3,I)*Quaternion(4,I))
        OMP_AMat_i(3,1) = 2.0D0*(Quaternion(1,I)*Quaternion(3,I)+Quaternion(2,I)*Quaternion(4,I))
        OMP_AMat_i(1,2) = 2.0D0*(Quaternion(1,I)*Quaternion(2,I)+Quaternion(3,I)*Quaternion(4,I))
        OMP_AMat_i(2,2) = -Quaternion(1,I)**2+Quaternion(2,I)**2-Quaternion(3,I)**2+Quaternion(4,I)**2
        OMP_AMat_i(3,2) = 2.0D0*(Quaternion(2,I)*Quaternion(3,I)-Quaternion(1,I)*Quaternion(4,I))
        OMP_AMat_i(1,3) = 2.0D0*(Quaternion(1,I)*Quaternion(3,I)-Quaternion(2,I)*Quaternion(4,I))
        OMP_AMat_i(2,3) = 2.0D0*(Quaternion(2,I)*Quaternion(3,I)+Quaternion(1,I)*Quaternion(4,I))
        OMP_AMat_i(3,3) = -Quaternion(1,I)**2-Quaternion(2,I)**2+Quaternion(3,I)**2+Quaternion(4,I)**2
        !  omega_i -> omega_b
        do K = 1,3
            OMP_Wb(K) = OMP_AMat_i(K,1)*W(1,I) + OMP_AMat_i(K,2)*W(2,I) + OMP_AMat_i(K,3)*W(3,I)
        end do
        !  omegadot_i -> omegadot_b
        do K = 1,3
            OMP_Wdotb(K) = OMP_AMat_i(K,1)*FM(1,I) + OMP_AMat_i(K,2)*FM(2,I) + OMP_AMat_i(K,3)*FM(3,I)
        end do
        !  Q, omega -> Qdot
        QuaternionDot(1) = 0.5D0 * ( OMP_Wb(1)*Quaternion(4,I) - OMP_Wb(2)*Quaternion(3,I) + OMP_Wb(3)*Quaternion(2,I))
        QuaternionDot(2) = 0.5D0 * ( OMP_Wb(1)*Quaternion(3,I) + OMP_Wb(2)*Quaternion(4,I) - OMP_Wb(3)*Quaternion(1,I))
        QuaternionDot(3) = 0.5D0 * (-OMP_Wb(1)*Quaternion(2,I) + OMP_Wb(2)*Quaternion(1,I) + OMP_Wb(3)*Quaternion(4,I))
        QuaternionDot(4) = 0.5D0 * (-OMP_Wb(1)*Quaternion(1,I) - OMP_Wb(2)*Quaternion(2,I) - OMP_Wb(3)*Quaternion(3,I))
        !  Q, omegadot, Qdot -> Qdotdot
        QuaternionDotDot(1) = ( OMP_Wdotb(1)*Quaternion(4,I) - OMP_Wdotb(2)*Quaternion(3,I) + OMP_Wdotb(3)*Quaternion(2,I))*0.5D0 + ( OMP_Wb(1)*QuaternionDot(4) -OMP_Wb(2)*QuaternionDot(3) + OMP_Wb(3)*QuaternionDot(2))*0.5D0
        QuaternionDotDot(2) = ( OMP_Wdotb(1)*Quaternion(3,I) + OMP_Wdotb(2)*Quaternion(4,I) - OMP_Wdotb(3)*Quaternion(1,I))*0.5D0 + ( OMP_Wb(1)*QuaternionDot(3) +OMP_Wb(2)*QuaternionDot(4) - OMP_Wb(3)*QuaternionDot(1))*0.5D0
        QuaternionDotDot(3) = (-OMP_Wdotb(1)*Quaternion(2,I) + OMP_Wdotb(2)*Quaternion(1,I) + OMP_Wdotb(3)*Quaternion(4,I))*0.5D0 + (-OMP_Wb(1)*QuaternionDot(2) +OMP_Wb(2)*QuaternionDot(1) + OMP_Wb(3)*QuaternionDot(4))*0.5D0
        QuaternionDotDot(4) = (-OMP_Wdotb(1)*Quaternion(1,I) - OMP_Wdotb(2)*Quaternion(2,I) - OMP_Wdotb(3)*Quaternion(3,I))*0.5D0 + (-OMP_Wb(1)*QuaternionDot(1) -OMP_Wb(2)*QuaternionDot(2) - OMP_Wb(3)*QuaternionDot(3))*0.5D0
        do K = 1,4
            Quaternion(K,I) = Quaternion(K,I) + QuaternionDot(K) * Dt +  QuaternionDotDot(K) * Dt * Dt /2.0D0
        end do
        norm = sqrt(Quaternion(1,I)**2 + Quaternion(2,I)**2 + Quaternion(3,I)**2 + Quaternion(4,I)**2)
        do K = 1,4
            Quaternion(K,I) = Quaternion(K,I)/norm
        end do
        
        !  refresh omega
        do K = 1,3
            W(K,I) = W(K,I) + FM(K,I) * Dt /2.0D0
        end do
    end do
    !$OMP END PARALLEL DO

    return
    end