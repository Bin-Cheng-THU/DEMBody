    !********************************************************************
    !     DEMBody 6.0
    !     ***********
    !
    !     N-body integrator of bond-assembly.
    !     -------------------------------
    !     @New method to solve the intgret: mid-step velocities verlet algoritm
    !     @Calculated x & xdot & w & Q from the previous values.
    !
    !********************************************************************
    subroutine intgrtAssembly()

    use global
    implicit none

    integer :: I,K,J
    real(8) :: Q(4),Qdot(4),Qdotdot(4)
    real(8) :: t_MatI(3,3),t_MatB(3,3)
    real(8) :: norm
    real(8) :: dist(3),distL
    integer :: index
    real(8) :: t_r(3)

    !  intgrt bond-assembly
    !$OMP PARALLEL DO PRIVATE(I,J,K,Q,Qdot,Qdotdot,norm,t_MatI,t_MatB)
    do I = 1,bondN
        !  refresh X & Xdot
        do K = 1,3
            bondMoveX(K,I) = bondXdot(K,I) * Dt + bondF(K,I) * Dt * Dt /2.0D0
            bondX(K,I) = bondX(K,I) + bondMoveX(K,I)
            bondXdot(K,I) = bondXdot(K,I) + bondF(K,I) * Dt /2.0D0
        end do

        !  refresh Q
        !  Q, omega -> Qdot
        Qdot(1) = ( bondWB(1,I)*bondQ(4,I) - bondWB(2,I)*bondQ(3,I) + bondWB(3,I)*bondQ(2,I))*0.5D0
        Qdot(2) = ( bondWB(1,I)*bondQ(3,I) + bondWB(2,I)*bondQ(4,I) - bondWB(3,I)*bondQ(1,I))*0.5D0
        Qdot(3) = (-bondWB(1,I)*bondQ(2,I) + bondWB(2,I)*bondQ(1,I) + bondWB(3,I)*bondQ(4,I))*0.5D0
        Qdot(4) = (-bondWB(1,I)*bondQ(1,I) - bondWB(2,I)*bondQ(2,I) - bondWB(3,I)*bondQ(3,I))*0.5D0
        !  Q, omegaDot, Qdot -> Qdotdot
        Qdotdot(1) = ( bondWdotB(1,I)*bondQ(4,I) - bondWdotB(2,I)*bondQ(3,I) + bondWdotB(3,I)*bondQ(2,I))*0.5D0 + ( bondWB(1,I)*Qdot(4) - bondWB(2,I)*Qdot(3) + bondWB(3,I)*Qdot(2))*0.5D0
        Qdotdot(2) = ( bondWdotB(1,I)*bondQ(3,I) + bondWdotB(2,I)*bondQ(4,I) - bondWdotB(3,I)*bondQ(1,I))*0.5D0 + ( bondWB(1,I)*Qdot(3) + bondWB(2,I)*Qdot(4) - bondWB(3,I)*Qdot(1))*0.5D0
        Qdotdot(3) = (-bondWdotB(1,I)*bondQ(2,I) + bondWdotB(2,I)*bondQ(1,I) + bondWdotB(3,I)*bondQ(4,I))*0.5D0 + (-bondWB(1,I)*Qdot(2) + bondWB(2,I)*Qdot(1) + bondWB(3,I)*Qdot(4))*0.5D0
        Qdotdot(4) = (-bondWdotB(1,I)*bondQ(1,I) - bondWdotB(2,I)*bondQ(2,I) - bondWdotB(3,I)*bondQ(3,I))*0.5D0 + (-bondWB(1,I)*Qdot(1) - bondWB(2,I)*Qdot(2) - bondWB(3,I)*Qdot(3))*0.5D0
        !  refresh Quaternion
        do K = 1,4
            Q(K) = bondQ(K,I) + Qdot(K) * Dt + Qdotdot(K) * Dt * Dt /2.0D0
        end do
        !  Normalize Quaternion
        norm = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2 + Q(4)**2)
        do K = 1,4
            Q(K) = Q(K)/norm
            bondQ(K,I) = Q(K)
        end do

        call attitudeQ2M(Q,t_MatI,t_MatB)

        !  refresh moveMat and matrix
        do K = 1,3
            do J = 1,3
                bondMoveMat(K,J,I) = t_MatB(K,J) - bondMatB(K,J,I)
                bondMatB(K,J,I) = t_MatB(K,J)
                bondMatI(K,J,I) = t_MatI(K,J)
            end do
        end do
        
        !  refresh Angular velocity in body frame
        do K = 1,3
            bondWB(K,I) = bondWB(K,I) + bondWdotB(K,I) * Dt /2.0D0
        end do
        !  refresh Angular velocity in inertia frame
        do K = 1,3
            bondW(K,I) = bondMatB(K,1,I)*bondWB(1,I) + bondMatB(K,2,I)*bondWB(2,I) + bondMatB(K,3,I)*bondWB(3,I)
        end do
    end do
    !$OMP END PARALLEL DO

    !  intgrt particles
    !$OMP PARALLEL DO PRIVATE(I,J,K,index,dist,distL,t_r)
    do I = 1,N
        !  refresh X & W
        index = bondTag(I)
        do K = 1,3
            dist(K) = bondMoveX(K,index) + bondMoveMat(K,1,index)*bondParticles(1,I) + bondMoveMat(K,2,index)*bondParticles(2,I) + bondMoveMat(K,3,index)*bondParticles(3,I)
            XT(K,I) = XT(K,I) + dist(K)            
            X(K,I) = X(K,I) + dist(K)
            W(K,I) = bondW(K,index)
        end do
#ifdef LatticeSearch        
        distL = XT(1,I)*XT(1,I) + XT(2,I)*XT(2,I) + XT(3,I)*XT(3,I)
        if (distL .GE. verlet) then
            refreshLattice = .true.
        end if
#endif LatticeSearch
        !  refresh R
        do K = 1,3
            t_r(K) = bondMatB(K,1,index)*bondParticles(1,I) + bondMatB(K,2,index)*bondParticles(2,I) + bondMatB(K,3,index)*bondParticles(3,I)
        end do
        !  refresh Xdot
        Xdot(1,I) = bondXdot(1,index) + W(2,I)*t_r(3) - W(3,I)*t_r(2)
        Xdot(2,I) = bondXdot(2,index) + W(3,I)*t_r(1) - W(1,I)*t_r(3)
        Xdot(3,I) = bondXdot(3,index) + W(1,I)*t_r(2) - W(2,I)*t_r(1)
    end do
    !$OMP END PARALLEL DO
    
    end