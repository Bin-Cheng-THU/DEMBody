    !********************************************************************
    !     DEMBody 7.0
    !     ***********
    !
    !     Generate mesh and lattice based on partition of biParticles.
    !     -------------------------------------------
    !     @Allocate particles into Paralle Lattice.
    !     @First, put the particle into inner region.
    !     @Second, put it into outer region, i.e., neighbor Parallel Lattice.
    !     @To speed up the DEM compution, we just use Top Left lattice to memory
    !     @particle IDs, and use Lower Right lattice to conduct force.
    !     @Tests show that half-search method is not suitable for this question,
    !     @so we change the code to full-search method!!!    
    !
    !     @Strategy: allocate particles into inner Lattices, so for particles in Lattice
    !     @     A, you just need to search contact of particles in Lattice A and Neighbor
    !     @     lattices (The neighbor lattices are stored in DEM)
    !
    !     @Allocate particles into Gravity Lattice.
    !     @Memory the total mass and mass center of particles in Gravity Lattices.
    !
    !********************************************************************
    subroutine latticeGenerateBiDisperseMix()

    use global
    use omp_lib
    implicit none

    integer :: I,J,K
    integer :: Tag
    logical :: Flag
    real(8) :: positionD(3),positionU(3)
    integer :: color
    real(8) :: o1,o2   
    type(biDisperseLattice),pointer :: Temp
    type(biDisperseLattice),pointer :: TempH  

    !o1 = omp_get_wtime()
    !$OMP PARALLEL DO PRIVATE(I)
    do I = 1,biDisperseNum
        !  calculate linklist employing X & LatDx
        MixLinklist(I) = floor((biDisperseX(1,I)+MixLatMx)/MixLatDx) + 1 + floor((biDisperseX(2,I)+MixLatMy)/MixLatDy)*MixLatNx+ floor((biDisperseX(3,I)+MixLatMz)/MixLatDz)*MixLatNx*MixLatNy
    end do
    !$OMP END PARALLEL DO
    !o2 = omp_get_wtime()
    !write(*,*) "MeshGenerate",(o2-o1)


    !o1 = omp_get_wtime()
    !  initialized Parallel Lattice
    !$OMP PARALLEL DO PRIVATE(I,Temp) SCHEDULE(GUIDED)
    do I = 1,MixLatNum
        !  reset Inner quantities    
        do while(associated(MixIDInner(I)%next))
            Temp => MixIDInner(I)%next
            MixIDInner(I)%next => MixIDInner(I)%next%next
            deallocate(Temp)
        end do
        MixIDInner(I)%No = 0
        !  assign tail of MixIDInner to tail
        MixtailInner(I)%next => MixIDInner(I)
    end do
    !$OMP END PARALLEL DO
    !o2 = omp_get_wtime()
    !write(*,*) "Lattice Empty",(o2-o1)
      
    !o1 = omp_get_wtime()
    ! $OMP PARALLEL DO PRIVATE(I,J,K,Tag,Flag) SCHEDULE(GUIDED)
    do I = 1,biDisperseNum
        
        !################         Part 1          ###################        
        Tag = MixLinklist(I)

        !****************
        !Note that in some cases the particles must be in the computational region, then
        !this check will be unnecessary, so we remove this part in Confined cases. 
        !****************
        !  insert into Lattice's inner region
        MixIDInner(Tag)%No = MixIDInner(Tag)%No + 1
        !TempH => MixIDInner(Tag)
        !do while(associated(TempH%next))
        !    TempH => TempH%next
        !end do
        allocate(Temp)
        Temp = biDisperseLattice(I,NULL())
        MixtailInner(Tag)%next%next => Temp
        MixtailInner(Tag)%next => Temp                
    end do
    ! $OMP END PARALLEL DO
    !o2 = omp_get_wtime()
    !write(*,*) "Lattice Assign",(o2-o1)
    
    end