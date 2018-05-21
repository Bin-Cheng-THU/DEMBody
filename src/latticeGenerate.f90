    !********************************************************************
    !     DEMBody 4.0
    !     ***********
    !
    !     Generate lattice based on partition of Particles.
    !     -------------------------------------------
    !     @Allocate particles into Paralle Lattice.
    !     @First, put the particle into inner region.
    !     @Second, put it into outer region, i.e., neighbor Parallel Lattice.
    !     @To speed up the DEM compution, we just use Top Left lattice to memory
    !     @particle IDs, and use Lower Right lattice to conduct force.
    !
    !     @Allocate particles into Gravity Lattice.
    !     @Memory the total mass and mass center of particles in Gravity Lattices.
    !
    !********************************************************************
    subroutine latticeGenerate()

    use global
    use omp_lib
    implicit none

    integer :: I,J,K
    integer :: Tag
    logical :: Flag
    real(8) :: positionD(3),positionU(3)
    integer :: color
    real(8) :: tmp_Mass,tmp_MassCenter(3)
    real(8) :: o1,o2
    
    !o1 = omp_get_wtime()
    !  initialized Parallel Lattice
    !$OMP PARALLEL DO PRIVATE(I)
    do I = 1,LatNum  
        !  reset Inner & Outer quantities
        DEM(I)%NoInner = 0               
        DEM(I)%NoOuter = 0               
        DEM(I)%IDInner = 0               
        DEM(I)%IDOuter = 0
    end do
    !$OMP END PARALLEL DO
    !o2 = omp_get_wtime()
    !write(*,*) (o2-o1)
       
    !o1 = omp_get_wtime()
    !$OMP PARALLEL DO PRIVATE(I,J,K,Tag,Flag,positionD,positionU) SCHEDULE(GUIDED)
    do I = 1,N
        
        !################         Part 1          ###################        
        Tag = floor((X(1,I)+Mx)/LatDx) + 1 + floor((X(2,I)+My)/LatDy)*LatNx + floor((X(3,I)+Mz)/LatDz)*(LatNx*LatNy)

        !****************
        !Note that if the particle is out of computational region, the Tag
        !may be beyond Normal Parallel Lattice, which would generat error when put into
        !DEM array; so we must check the position of the particle.
        !****************

        if (Tag.GE.1 .AND. Tag.LE.LatNum) then
            
            Flag = .True.

            !  check position of the particle
            positionD = DEM(Tag)%PositionD
            positionU = DEM(Tag)%PositionU

            !  out of the region
            if (X(1,I).LT.positionD(1) .OR. X(1,I).GT.positionU(1) .OR. X(2,I).LT.positionD(2) .OR. X(2,I).GT.positionU(2) &
            .OR. X(3,I).LT.positionD(3) .OR. X(3,I).GT.positionU(3)) then
                Flag = .False.
            end if

            if (Flag) then
                !  insert into Lattice's inner region
                DEM(Tag)%NoInner = DEM(Tag)%NoInner + 1
                DEM(Tag)%IDInner(DEM(Tag)%NoInner) = I

                !  insert into Lattice's outer region
                !  Left lattice
                if (DEM(Tag)%ID(1) .NE. 1) then
                    if (X(1,I) .LT. (positionD(1)+Dx)) then
                        DEM(Tag-1)%NoOuter = DEM(Tag-1)%NoOuter + 1
                        DEM(Tag-1)%IDOuter(DEM(Tag-1)%NoOuter) = I
                    end if
                end if

                !  Left Upper lattice
                if (DEM(Tag)%ID(1).NE.1 .AND. DEM(Tag)%ID(2).NE.LatNy) then
                    if (X(1,I).LT.(positionD(1)+Dx) .AND. X(2,I).GT.(positionU(2)-Dy)) then
                        DEM(Tag-1+LatNx)%NoOuter = DEM(Tag-1+LatNx)%NoOuter + 1
                        DEM(Tag-1+LatNx)%IDOuter(DEM(Tag-1+LatNx)%NoOuter) = I
                    end if
                end if

                !  Upper lattice
                if (DEM(Tag)%ID(2) .NE. LatNy) then
                    if (X(2,I) .GT. (positionU(2)-Dy)) then
                        DEM(Tag+LatNx)%NoOuter = DEM(Tag+LatNx)%NoOuter + 1
                        DEM(Tag+LatNx)%IDOuter(DEM(Tag+LatNx)%NoOuter) = I
                    end if
                end if            

                !  Right Upper lattice
                if (DEM(Tag)%ID(1).NE.LatNx .AND. DEM(Tag)%ID(2).NE.LatNy) then
                    if (X(1,I).GT.(positionU(1)-Dx) .AND. X(2,I).GT.(positionU(2)-Dy)) then
                        DEM(Tag+1+LatNx)%NoOuter = DEM(Tag+1+LatNx)%NoOuter + 1
                        DEM(Tag+1+LatNx)%IDOuter(DEM(Tag+1+LatNx)%NoOuter) = I
                    end if
                end if

                !  Center Covered lattice
                if (DEM(Tag)%ID(3) .NE. LatNz) then
                    if (X(3,I) .GT. (positionU(3)-Dz)) then
                        DEM(Tag+LatNx*LatNy)%NoOuter = DEM(Tag+LatNx*LatNy)%NoOuter + 1
                        DEM(Tag+LatNx*LatNy)%IDOuter(DEM(Tag+LatNx*LatNy)%NoOuter) = I
                    end if
                end if

                !  Left Covered lattice
                if (DEM(Tag)%ID(1).NE.1 .AND. DEM(Tag)%ID(3).NE.LatNz) then
                    if (X(1,I).LT.(positionD(1)+Dx) .AND. X(3,I).GT.(positionU(3)-Dz)) then
                        DEM(Tag-1+LatNx*LatNy)%NoOuter = DEM(Tag-1+LatNx*LatNy)%NoOuter + 1
                        DEM(Tag-1+LatNx*LatNy)%IDOuter(DEM(Tag-1+LatNx*LatNy)%NoOuter) = I
                    end if
                end if

                !  Left Upper Covered lattice
                if (DEM(Tag)%ID(1).NE.1 .AND. DEM(Tag)%ID(2).NE.LatNy .AND. DEM(Tag)%ID(3).NE.LatNz) then
                    if (X(1,I).LT.(positionD(1)+Dx) .AND. X(2,I).GT.(positionU(2)-Dy) .AND. X(3,I).GT.(positionU(3)-Dz)) then
                        DEM(Tag-1+LatNx+LatNx*LatNy)%NoOuter = DEM(Tag-1+LatNx+LatNx*LatNy)%NoOuter + 1
                        DEM(Tag-1+LatNx+LatNx*LatNy)%IDOuter(DEM(Tag-1+LatNx+LatNx*LatNy)%NoOuter) = I
                    end if
                end if

                !  Upper Covered lattice
                if (DEM(Tag)%ID(2).NE.LatNy .AND. DEM(Tag)%ID(3).NE.LatNz) then
                    if (X(2,I).GT.(positionU(2)-Dy) .AND. X(3,I).GT.(positionU(3)-Dz)) then
                        DEM(Tag+LatNx+LatNx*LatNy)%NoOuter = DEM(Tag+LatNx+LatNx*LatNy)%NoOuter + 1
                        DEM(Tag+LatNx+LatNx*LatNy)%IDOuter(DEM(Tag+LatNx+LatNx*LatNy)%NoOuter) = I
                    end if
                end if             

                !  Right Upper Covered lattice
                if (DEM(Tag)%ID(1).NE.LatNx .AND. DEM(Tag)%ID(2).NE.LatNy .AND. DEM(Tag)%ID(3).NE.LatNz) then
                    if (X(1,I).GT.(positionU(1)-Dx) .AND. X(2,I).GT.(positionU(2)-Dy) .AND. X(3,I).GT.(positionU(3)-Dz)) then
                        DEM(Tag+1+LatNx+LatNx*LatNy)%NoOuter = DEM(Tag+1+LatNx+LatNx*LatNy)%NoOuter + 1
                        DEM(Tag+1+LatNx+LatNx*LatNy)%IDOuter(DEM(Tag+1+LatNx+LatNx*LatNy)%NoOuter) = I
                    end if
                end if

                !  Right Covered lattice
                if (DEM(Tag)%ID(1).NE.LatNx .AND. DEM(Tag)%ID(3).NE.LatNz) then
                    if (X(1,I).GT.(positionU(1)-Dx) .AND. X(3,I).GT.(positionU(3)-Dz)) then
                        DEM(Tag+1+LatNx*LatNy)%NoOuter = DEM(Tag+1+LatNx*LatNy)%NoOuter + 1
                        DEM(Tag+1+LatNx*LatNy)%IDOuter(DEM(Tag+1+LatNx*LatNy)%NoOuter) = I
                    end if
                end if

                !  Right Below Covered lattice
                if (DEM(Tag)%ID(1).NE.LatNx .AND. DEM(Tag)%ID(2).NE.1 .AND. DEM(Tag)%ID(3).NE.LatNz) then
                    if (X(1,I).GT.(positionU(1)-Dx) .AND. X(2,I).LT.(positionD(2)+Dy) .AND. X(3,I).GT.(positionU(3)-Dz)) then
                        DEM(Tag+1-LatNx+LatNx*LatNy)%NoOuter = DEM(Tag+1-LatNx+LatNx*LatNy)%NoOuter + 1
                        DEM(Tag+1-LatNx+LatNx*LatNy)%IDOuter(DEM(Tag+1-LatNx+LatNx*LatNy)%NoOuter) = I
                    end if
                end if

                !  Below Covered lattice
                if (DEM(Tag)%ID(2).NE.1 .AND. DEM(Tag)%ID(3).NE.LatNz) then
                    if (X(2,I).LT.(positionD(2)+Dy) .AND. X(3,I).GT.(positionU(3)-Dz)) then
                        DEM(Tag-LatNx+LatNx*LatNy)%NoOuter = DEM(Tag-LatNx+LatNx*LatNy)%NoOuter + 1
                        DEM(Tag-LatNx+LatNx*LatNy)%IDOuter(DEM(Tag-LatNx+LatNx*LatNy)%NoOuter) = I
                    end if
                end if            

                !  Left Below Covered lattice
                if (DEM(Tag)%ID(1).NE.1 .AND. DEM(Tag)%ID(2).NE.1 .AND. DEM(Tag)%ID(3).NE.LatNz) then
                    if (X(1,I).LT.(positionD(1)+Dx) .AND. X(2,I).LT.(positionD(2)+Dy) .AND. X(3,I).GT.(positionU(3)-Dz)) then
                        DEM(Tag-1-LatNx+LatNx*LatNy)%NoOuter = DEM(Tag-1-LatNx+LatNx*LatNy)%NoOuter + 1
                        DEM(Tag-1-LatNx+LatNx*LatNy)%IDOuter(DEM(Tag-1-LatNx+LatNx*LatNy)%NoOuter) = I
                    end if
                end if
            end if
        end if
    end do
    !$OMP END PARALLEL DO
    !o2 = omp_get_wtime()
    !write(*,*) (o2-o1)
    
!!#ifdef self_gravity
!    !o1 = omp_get_wtime()
!    !$OMP PARALLEL DO PRIVATE(I,J,K,tmp_Mass,tmp_MassCenter)
!    !################         Part 2          ###################
!    do I = 1,LatNum
!        if (DEM(I)%NoInner .GE. 1) then
!            tmp_Mass = 0.0D0
!            tmp_MassCenter = 0.0D0
!            
!            do J = 1,DEM(I)%NoInner
!                do K = 1,3
!                    tmp_MassCenter(K) = tmp_MassCenter(K) + Body(DEM(I)%IDInner(J))*X(K,DEM(I)%IDInner(J))
!                end do
!                tmp_Mass = tmp_Mass + Body(DEM(I)%IDInner(J))
!            end do
!            
!            do K = 1,3
!                DEM(I)%MassCenter = tmp_MassCenter/tmp_Mass
!            end do
!            DEM(I)%Mass = tmp_Mass
!        else
!            DEM(I)%MassCenter = 0.0D0
!            DEM(I)%Mass = 0.0D0
!        end if
!    end do
!    !$OMP END PARALLEL DO
!    !o2 = omp_get_wtime()
!    !write(*,*) (o2-o1)
!    
!    !o1 = omp_get_wtime()
!    !$OMP PARALLEL DO PRIVATE(I,J,K,tmp_Mass,tmp_MassCenter)
!    !################         Part 3          ###################
!    do I = 1,GravNum
!        tmp_Mass = 0.0D0
!        tmp_MassCenter = 0.0D0
!
!        do J = 1,Gravity(I)%num
!            do K = 1,3
!                tmp_MassCenter(K) = tmp_MassCenter(K) + DEM(Gravity(I)%ID(J))%Mass*DEM(Gravity(I)%ID(J))%MassCenter(K)
!            end do
!            tmp_Mass = tmp_Mass + DEM(Gravity(I)%ID(J))%Mass
!        end do
!        
!        if (tmp_Mass .GT. 1e-10) then
!            Gravity(I)%MassCenter = tmp_MassCenter/tmp_Mass
!            Gravity(I)%Mass = tmp_Mass
!        else
!            Gravity(I)%MassCenter = 0.0D0
!            Gravity(I)%Mass = 0.0D0
!        end if
!    end do
!    !$OMP END PARALLEL DO
!    !o2 = omp_get_wtime()
!    !write(*,*) (o2-o1)
!!#endif    

    !****************
    !Test Parallel Lattice Method using color visualization.
    !****************
    
    !do I = 1,GravNum
    !    do J = 1,Gravity(I)%num
    !        do K = 1,DEM(Gravity(I)%ID(J))%NoInner
    !            ParallelLatticeColor(1,DEM(Gravity(I)%ID(J))%IDInner(K)) = I
    !        end do
    !    end do
    !end do
    !                
    !open(100,FILE='../Data/DEM.txt')
    !write(100,*) 'X',',','Y',',','Z',',','R',',','Color'
    !do I = 1,N
    !    write(100,'(F15.5,A2,F15.5,A2,F15.5,A2,F15.5,A2,I5)') X(1,I),',',X(2,I),',',X(3,I),',',R(I),',',ParallelLatticeColor(1,I)
    !end do
    !
    !stop

    end