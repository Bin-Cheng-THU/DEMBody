    !********************************************************************
    !     DEMBody 5.0
    !     ***********
    !
    !     Generate lattice based on partition of Particles.
    !     -------------------------------------------
    !     @Allocate particles into Paralle Lattice.
    !     @First, put the particle into inner region.
    !     @Second, put it into outer region, i.e., neighbor Parallel Lattice.
    !     @To speed up the DEM compution, we just use Top Left lattice to memory
    !     @particle IDs, and use Lower Right lattice to conduct force.
    !     @Tests show that half-search method is not suitable for this question,
    !     @so we change the code to full-search method.    
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
    real(8) :: o1,o2   
    type(Neighbor),pointer :: Temp
    type(Neighbor),pointer :: TempH  
    
    !o1 = omp_get_wtime()
    !  initialized Parallel Lattice
    !$OMP PARALLEL DO PRIVATE(I,Temp) SCHEDULE(GUIDED)
    do I = 1,LatNum
        !  reset Inner quantities    
        do while(associated(IDInner(I)%next))
            Temp => IDInner(I)%next
            IDInner(I)%next => IDInner(I)%next%next
            deallocate(Temp)
        end do
        IDInner(I)%No = 0
        !  assign tail of Inner to tail
        tailInner(I)%next => IDInner(I)
    end do
    !$OMP END PARALLEL DO
    !o2 = omp_get_wtime()
    !write(*,*) "Lattice Empty",(o2-o1)
       
    !o1 = omp_get_wtime()
    ! $OMP PARALLEL DO PRIVATE(I,J,K,Tag,Flag,positionD,positionU) SCHEDULE(GUIDED)
    do I = 1,N
        
        !################         Part 1          ###################        
        !Tag = floor((X(1,I)+LatMx)/LatDx) + 1 + floor((X(2,I)+LatMy)/LatDy)*LatNx + floor((X(3,I)+LatMz)/LatDz)*(LatNx*LatNy)
        Tag = Linklist(I)

#ifdef nConfined
        !****************
        !Note that if the particle is out of computational region, the Tag
        !may be beyond Normal Parallel Lattice, which would generat error when put into
        !DEM array; so we must check the position of the particle.
        !****************
        if (Tag.GE.1 .AND. Tag.LE.LatNum) then
            
            Flag = .true.

            !  check position of the particle
            positionD = DEM(Tag)%PositionD
            positionU = DEM(Tag)%PositionU

            !  out of the region
            if (X(1,I).LT.positionD(1) .OR. X(1,I).GT.positionU(1) .OR. X(2,I).LT.positionD(2) .OR. X(2,I).GT.positionU(2) &
            .OR. X(3,I).LT.positionD(3) .OR. X(3,I).GT.positionU(3)) then
                Flag = .false.
            end if

            if (Flag) then
                !  insert into Lattice's inner region
                IDInner(Tag)%No = IDInner(Tag)%No + 1
                !TempH => IDInner(Tag)
                !do while(associated(TempH%next))
                !    TempH => TempH%next
                !end do        
                allocate(Temp)
                Temp = Neighbor(I,NULL())
                tailInner(Tag)%next%next => Temp
                tailInner(Tag)%next => Temp    
            end if         
        end if
#elif Confined
        !****************
        !Note that in some cases the particles must be in the computational region, then
        !this check will be unnecessary, so we remove this part in Confined cases. 
        !****************

        !  insert into Lattice's inner region
        IDInner(Tag)%No = IDInner(Tag)%No + 1
        !TempH => IDInner(Tag)
        !do while(associated(TempH%next))
        !    TempH => TempH%next
        !end do        
        allocate(Temp)
        Temp = Neighbor(I,NULL())
        tailInner(Tag)%next%next => Temp
        tailInner(Tag)%next => Temp         
#endif        
    end do
    ! $OMP END PARALLEL DO
    !o2 = omp_get_wtime()
    !write(*,*) "Lattice Assign",(o2-o1)
    
    !  reset Position array
    !$OMP PARALLEL DO PRIVATE(I,K)
    do I = 1,N
        do K = 1,3
            XT(K,I) = 0.0D0
        end do
    end do
    !$OMP END PARALLEL DO
    
    !  reset Lattice flag
    refreshNum = refreshNum + 1
    refreshLattice = .false.

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
    !close(100)
    !
    !stop
    
    !****************
    !Test Parallel Lattice Method using files.
    !****************
!#ifdef ArrayStore
!    open(100,FILE="Inner.txt")
!    do J = 1,LatNum
!        write(100,'(I5)',advance='no') DEM(J)%NoInner
!        do I = 1,DEM(J)%NoInner
!            write(100,'(I5)',advance='no') DEM(J)%IDInner(I)
!        end do
!        write(100,*)
!    end do
!    close(100)
!    open(100,FILE="Outer.txt")
!    do J = 1,LatNum
!        write(100,'(I5)',advance='no') DEM(J)%NoOuter
!        do I = 1,DEM(J)%NoOuter
!            write(100,'(I5)',advance='no') DEM(J)%IDOuter(I)
!        end do
!        write(100,*)
!    end do
!    close(100)
!#elif LinklistStore
!    open(100,FILE="Inner.txt")
!    do J = 1,LatNum
!        TempH => IDInner(J)
!        write(100,'(I5)',advance='no') TempH%No
!        do while(associated(TempH%next))
!            TempH => TempH%next
!            write(100,'(I5)',advance='no') TempH%No
!        end do 
!        write(100,*)
!    end do
!    close(100)
!    open(100,FILE="Outer.txt")
!    do J = 1,LatNum
!        TempH => IDOuter(J)
!        write(100,'(I5)',advance='no') TempH%No
!        do while(associated(TempH%next))
!            TempH => TempH%next
!            write(100,'(I5)',advance='no') TempH%No
!        end do 
!        write(100,*)
!    end do
!    close(100)
!#endif
!
!    pause
    
    end