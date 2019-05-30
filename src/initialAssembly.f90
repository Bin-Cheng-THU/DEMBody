    !********************************************************************
    !     DEMBody 6.0
    !     ***********
    !
    !     Conduct particle-kinematics using assembly-kinematics.
    !     -------------------------------
    !     
    !     @Read data for assmbly informations
    !     @Refresh angular velocity, position and velocity
    !
    !********************************************************************
    subroutine initialAssembly()
    
    use global
    implicit none
    
    integer :: I,J,K
    integer :: index
    real(8) :: t_r(3)
    real(8) :: t_quaternion(4),t_MatI(3,3),t_MatB(3,3)
    integer :: count
    
    allocate(bondNumN(bondN))
    allocate(bondX(3,bondN))
    allocate(bondXdot(3,bondN))
    allocate(bondW(3,bondN))
    allocate(bondBody(bondN))
    allocate(bondInertia(3,bondN))
    allocate(bondF(3,bondN))
    allocate(bondFM(3,bondN))
    allocate(bondQ(4,bondN))
    allocate(bondMatI(3,3,bondN))
    allocate(bondMatB(3,3,bondN))
    allocate(bondWB(3,bondN))
    allocate(bondWdotB(3,bondN))
    allocate(bondMoveX(3,bondN))
    allocate(bondMoveMat(3,3,bondN))
    allocate(bondRange(2,bondN))
    
    count = 1
    do I = 1,bondN
        read (3000,*) bondNumN(I),bondBody(I),(bondInertia(K,I),K=1,3),(bondX(K,I),K=1,3),(bondXdot(K,I),K=1,3),(bondW(K,I),K=1,3),(bondQ(K,I),K=1,4)
        !  Tag & Range
        bondRange(1,I) = count
        do J = 1,bondNumN(I)
            bondTag(count) = I
            count = count + 1
        end do
        bondRange(2,I) = count - 1
        !  Q & Matrix
        do K = 1,4
            t_quaternion(K) = bondQ(K,I)
        end do
        call attitudeQ2M(t_quaternion,t_MatI,t_MatB)
        do K = 1,3
            do J = 1,3
                bondMatI(K,J,I) = t_MatI(K,J)
                bondMatB(K,J,I) = t_MatB(K,J)
            end do
        end do
        !  Wb
        do K = 1,3
            bondWB(K,I) = bondMatI(K,1,I)*bondW(1,I) + bondMatI(K,2,I)*bondW(2,I) + bondMatI(K,3,I)*bondW(3,I)
        end do
    end do
    
    !$OMP PARALLEL DO PRIVATE(I,K,index,t_r)
    do I = 1,N
        index = bondTag(I)
        !  refresh W
        do K = 1,3
            W(K,I) = bondW(K,index)
        end do
        !  refresh X
        do K = 1,3
            t_r(K) = bondMatB(K,1,index)*bondParticles(1,I) + bondMatB(K,2,index)*bondParticles(2,I) + bondMatB(K,3,index)*bondParticles(3,I)
        end do
        do K = 1,3
            X(K,I) = bondX(K,index) + t_r(K)
        end do
        !  refresh Xdot
        Xdot(1,I) = bondXdot(1,index) + W(2,I)*t_r(3) - W(3,I)*t_r(2)
        Xdot(2,I) = bondXdot(2,index) + W(3,I)*t_r(1) - W(1,I)*t_r(3)
        Xdot(3,I) = bondXdot(3,index) + W(1,I)*t_r(2) - W(2,I)*t_r(1)
    end do
    !$OMP END PARALLEL DO
    
    end subroutine