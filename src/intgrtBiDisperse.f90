    !********************************************************************
    !     DEMBody 7.0
    !     ***********
    !
    !     N-body integrator of BiDisperse.
    !     -------------------------------
    !     @Refresh the position of BiDisperse in T
    !     @Refresh the velocity/angular velocity of BiDisperse in T/2
    !     @Refresh the quaternion of BiDisperse in T
    !
    !     @Note that we need to calculate the displacement of particles to check whether
    !     @we need to refresh the Lattice for biDisperse particles
    !
    !********************************************************************
    subroutine intgrtBiDisperse()

    use global
    implicit none

    integer :: I,K
    real(8) :: dist(3),distL
    real(8) :: biDisperseQdot(4),biDisperseQdotdot(4)
    real(8) :: norm
    real(8) :: biDisperseAMat_i(3,3)
    real(8) :: biDisperseWb(3), biDisperseWdotb(3)

    !  refresh X, V
    !$OMP PARALLEL DO PRIVATE(I,K,dist,distL)
    do I = 1,biDisperseNum
        do K = 1,3
            dist(K) = biDisperseXdot(K,I) * Dt +  biDisperseF(K,I) * Dt * Dt /2.0D0
            biDisperseXT(K,I) = biDisperseXT(K,I) + dist(K)            
            biDisperseX(K,I) = biDisperseX(K,I) + dist(K)
            biDisperseXdot(K,I) = biDisperseXdot(K,I) + biDisperseF(K,I) * Dt /2.0D0
        end do
#ifdef LatticeSearch
        distL = biDisperseXT(1,I)*biDisperseXT(1,I) + biDisperseXT(2,I)*biDisperseXT(2,I) + biDisperseXT(3,I)*biDisperseXT(3,I)
        if (distL .GE. verlet) then
            refreshLattice = .true.
        end if
#endif LatticeSearch
    end do
    !$OMP END PARALLEL DO
    
    if (isQuaternion) then
        !$OMP PARALLEL DO PRIVATE(I,K,biDisperseAMat_i,biDisperseWb,biDisperseWdotb,biDisperseQdot,biDisperseQdotdot,norm)
        do I = 1,biDisperseNum
            !  refresh Q
            !  Q -> AMat_i
            biDisperseAMat_i(1,1) = biDisperseQ(1,I)**2-biDisperseQ(2,I)**2-biDisperseQ(3,I)**2+biDisperseQ(4,I)**2
            biDisperseAMat_i(2,1) = 2.0D0*(biDisperseQ(1,I)*biDisperseQ(2,I)-biDisperseQ(3,I)*biDisperseQ(4,I))
            biDisperseAMat_i(3,1) = 2.0D0*(biDisperseQ(1,I)*biDisperseQ(3,I)+biDisperseQ(2,I)*biDisperseQ(4,I))
            biDisperseAMat_i(1,2) = 2.0D0*(biDisperseQ(1,I)*biDisperseQ(2,I)+biDisperseQ(3,I)*biDisperseQ(4,I))
            biDisperseAMat_i(2,2) = -biDisperseQ(1,I)**2+biDisperseQ(2,I)**2-biDisperseQ(3,I)**2+biDisperseQ(4,I)**2
            biDisperseAMat_i(3,2) = 2.0D0*(biDisperseQ(2,I)*biDisperseQ(3,I)-biDisperseQ(1,I)*biDisperseQ(4,I))
            biDisperseAMat_i(1,3) = 2.0D0*(biDisperseQ(1,I)*biDisperseQ(3,I)-biDisperseQ(2,I)*biDisperseQ(4,I))
            biDisperseAMat_i(2,3) = 2.0D0*(biDisperseQ(2,I)*biDisperseQ(3,I)+biDisperseQ(1,I)*biDisperseQ(4,I))
            biDisperseAMat_i(3,3) = -biDisperseQ(1,I)**2-biDisperseQ(2,I)**2+biDisperseQ(3,I)**2+biDisperseQ(4,I)**2
            !  omega_i -> omega_b
            do K = 1,3
                biDisperseWb(K) = biDisperseAMat_i(K,1)*biDisperseW(1,I) + biDisperseAMat_i(K,2)*biDisperseW(2,I) + biDisperseAMat_i(K,3)*biDisperseW(3,I)
            end do
            !  omegadot_i -> omegadot_b
            do K = 1,3
                biDisperseWdotb(K) = biDisperseAMat_i(K,1)*biDisperseFM(1,I) + biDisperseAMat_i(K,2)*biDisperseFM(2,I) + biDisperseAMat_i(K,3)*biDisperseFM(3,I)
            end do
            !  Q, omega -> Qdot
            biDisperseQdot(1) = 0.5D0 * ( biDisperseWb(1)*biDisperseQ(4,I) - biDisperseWb(2)*biDisperseQ(3,I) + biDisperseWb(3)*biDisperseQ(2,I))
            biDisperseQdot(2) = 0.5D0 * ( biDisperseWb(1)*biDisperseQ(3,I) + biDisperseWb(2)*biDisperseQ(4,I) - biDisperseWb(3)*biDisperseQ(1,I))
            biDisperseQdot(3) = 0.5D0 * (-biDisperseWb(1)*biDisperseQ(2,I) + biDisperseWb(2)*biDisperseQ(1,I) + biDisperseWb(3)*biDisperseQ(4,I))
            biDisperseQdot(4) = 0.5D0 * (-biDisperseWb(1)*biDisperseQ(1,I) - biDisperseWb(2)*biDisperseQ(2,I) - biDisperseWb(3)*biDisperseQ(3,I))
            !  Q, omegadot, Qdot -> Qdotdot
            biDisperseQdotdot(1) = ( biDisperseWdotb(1)*biDisperseQ(4,I) - biDisperseWdotb(2)*biDisperseQ(3,I) + biDisperseWdotb(3)*biDisperseQ(2,I))*0.5D0 + ( biDisperseWb(1)*biDisperseQdot(4) - biDisperseWb(2)*biDisperseQdot(3) + biDisperseWb(3)*biDisperseQdot(2))*0.5D0
            biDisperseQdotdot(2) = ( biDisperseWdotb(1)*biDisperseQ(3,I) + biDisperseWdotb(2)*biDisperseQ(4,I) - biDisperseWdotb(3)*biDisperseQ(1,I))*0.5D0 + ( biDisperseWb(1)*biDisperseQdot(3) + biDisperseWb(2)*biDisperseQdot(4) - biDisperseWb(3)*biDisperseQdot(1))*0.5D0
            biDisperseQdotdot(3) = (-biDisperseWdotb(1)*biDisperseQ(2,I) + biDisperseWdotb(2)*biDisperseQ(1,I) + biDisperseWdotb(3)*biDisperseQ(4,I))*0.5D0 + (-biDisperseWb(1)*biDisperseQdot(2) + biDisperseWb(2)*biDisperseQdot(1) + biDisperseWb(3)*biDisperseQdot(4))*0.5D0
            biDisperseQdotdot(4) = (-biDisperseWdotb(1)*biDisperseQ(1,I) - biDisperseWdotb(2)*biDisperseQ(2,I) - biDisperseWdotb(3)*biDisperseQ(3,I))*0.5D0 + (-biDisperseWb(1)*biDisperseQdot(1) - biDisperseWb(2)*biDisperseQdot(2) - biDisperseWb(3)*biDisperseQdot(3))*0.5D0
            do K = 1,4
                biDisperseQ(K,I) = biDisperseQ(K,I) + biDisperseQdot(K) * Dt + biDisperseQdotdot(K) * Dt * Dt /2.0D0
            end do
            norm = sqrt(biDisperseQ(1,I)**2 + biDisperseQ(2,I)**2 + biDisperseQ(3,I)**2 + biDisperseQ(4,I)**2)
            do K = 1,4
                biDisperseQ(K,I) = biDisperseQ(K,I)/norm
            end do
        end do
        !$OMP END PARALLEL DO
    end if

    !$OMP PARALLEL DO PRIVATE(I,K)
    do I = 1,biDisperseNum
        !  refresh omega
        do K = 1,3
            biDisperseW(K,I) = biDisperseW(K,I) + biDisperseFM(K,I) * Dt /2.0D0
        end do
    end do
    !$OMP END PARALLEL DO
    
    return
    end