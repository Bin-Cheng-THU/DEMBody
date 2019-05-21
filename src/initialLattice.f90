    !********************************************************************
    !     DEMBody 5.2
    !     ***********
    !
    !     Initialization of DEM Lattice.
    !     ---------------------------------
    !
    !********************************************************************  
    SUBROUTINE initialLattice()

    use global
    implicit none
    
    integer I,J,K,L
    integer :: idx,idy,idz
    integer :: idgx,idgy,idgz
    real(8) :: lx,ly,lz
    
    integer :: idLattice
    type(PeriodicLattice),pointer :: Temp
    type(PeriodicLattice),pointer :: Tail
    real(8) :: ostart,oend
    
    real(8),allocatable :: test_lattice(:,:)
    integer,allocatable :: test_color(:)
    integer :: count
    integer :: test_index

#ifdef LatticeSearch
    !  Initialize Parallel Lattice.
    write(*,*) "< Parallel Lattice initializing..."
    
    LatNum = LatNx*LatNy*LatNz             !  Paraller Number
    allocate(DEM(LatNum))
    do I = 1,LatNum
        idz = int((I-1)/(LatNx*LatNy))+1
        idy = int((I-(idz-1)*(LatNx*LatNy)-1)/LatNx)+1
        idx = I-(idz-1)*(LatNx*LatNy)-(idy-1)*LatNx
        lx = dble(idx-1)*LatDx - LatMx
        ly = dble(idy-1)*LatDy - LatMy
        lz = dble(idz-1)*LatDz - LatMz
        !DEM(I)%ID(1) = idx                 !  IDx of Lattice
        !DEM(I)%ID(2) = idy                 !  IDy of Lattice
        !DEM(I)%ID(3) = idz                 !  IDz of Lattice
#ifdef nConfined        
        DEM(I)%PositionD(1) = lx           !  PositionDx of Lattice
        DEM(I)%PositionD(2) = ly           !  PositionDy of Lattice
        DEM(I)%PositionD(3) = lz           !  PositionDz of Lattice
        DEM(I)%PositionU(1) = lx + LatDx   !  PositionUx of Lattice
        DEM(I)%PositionU(2) = ly + LatDy   !  PositionUy of Lattice
        DEM(I)%PositionU(3) = lz + LatDz   !  PositionUz of Lattice
#endif      
        DEM(I)%NeighborID = 0
        !  Left lattice
        if (idx .NE. 1) then
            DEM(I)%NeighborID(1) = I - 1
        end if
        !  Left Upper lattice
        if (idx.NE.1 .AND. idy.NE.LatNy) then
            DEM(I)%NeighborID(2) = I - 1 + LatNx
        end if
        !  Upper lattice
        if (idy .NE. LatNy) then
            DEM(I)%NeighborID(3) = I + LatNx
        end if
        !  Right Upper lattice
        if (idx.NE.LatNx .AND. idy.NE.LatNy) then
            DEM(I)%NeighborID(4) = I + 1 + LatNx
        end if
        !  Center Covered lattice
        if (idz .NE. LatNz) then
            DEM(I)%NeighborID(5) = I + LatNx*LatNy
        end if
        !  Left Covered lattice
        if (idx.NE.1 .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(6) = I - 1 + LatNx*LatNy
        end if
        !  Left Upper Covered lattice
        if (idx.NE.1 .AND. idy.NE.LatNy .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(7) = I - 1 + LatNx + LatNx*LatNy
        end if
        !  Upper Covered lattice
        if (idy.NE.LatNy .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(8) = I + LatNx + LatNx*LatNy
        end if
        !  Right Upper Covered lattice
        if (idx.NE.LatNx .AND. idy.NE.LatNy .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(9) = I + 1 + LatNx + LatNx*LatNy
        end if
        !  Right Covered lattice
        if (idx.NE.LatNx .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(10) = I + 1 + LatNx*LatNy
        end if
        !  Right Below Covered lattice
        if (idx.NE.LatNx .AND. idy.NE.1 .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(11) = I + 1 - LatNx + LatNx*LatNy
        end if
        !  Below Covered lattice
        if (idy.NE.1 .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(12) = I - LatNx + LatNx*LatNy
        end if
        !  Left Below Covered lattice
        if (idx.NE.1 .AND. idy.NE.1 .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(13) = I - 1 - LatNx + LatNx*LatNy
        end if
        !  Right lattice
        if (idx .NE. LatNx) then
            DEM(I)%NeighborID(14) = I + 1
        end if
        !  Right Below lattice
        if (idx.NE.LatNx .AND. idy.NE.1) then
            DEM(I)%NeighborID(15) = I + 1 - LatNx
        end if
        !  Below lattice
        if (idy .NE. 1) then
            DEM(I)%NeighborID(16) = I - LatNx
        end if
        !  Left Below lattice
        if (idx.NE.1 .AND. idy.NE.1) then
            DEM(I)%NeighborID(17) = I - 1 - LatNx
        end if
        !  Center Under lattice
        if (idz .NE. 1) then
            DEM(I)%NeighborID(18) = I - LatNx*LatNy
        end if
        !  Right Under lattice
        if (idx.NE.LatNx .AND. idz.NE.1) then
            DEM(I)%NeighborID(19) = I + 1 - LatNx*LatNy
        end if
        !  Right Below Under lattice
        if (idx.NE.LatNx .AND. idy.NE.1 .AND. idz.NE.1) then
            DEM(I)%NeighborID(20) = I + 1 - LatNx - LatNx*LatNy
        end if
        !  Below Under lattice
        if (idy.NE.1 .AND. idz.NE.1) then
            DEM(I)%NeighborID(21) = I - LatNx - LatNx*LatNy
        end if
        !  Left Below Under lattice
        if (idx.NE.1 .AND. idy.NE.1 .AND. idz.NE.1) then
            DEM(I)%NeighborID(22) = I - 1 - LatNx - LatNx*LatNy
        end if
        !  Left Under lattice
        if (idx.NE.1 .AND. idz.NE.1) then
            DEM(I)%NeighborID(23) = I - 1 - LatNx*LatNy
        end if
        !  Left Upper Under lattice
        if (idx.NE.1 .AND. idy.NE.LatNy .AND. idz.NE.1) then
            DEM(I)%NeighborID(24) = I - 1 + LatNx - LatNx*LatNy
        end if
        !  Upper Under lattice
        if (idy.NE.LatNy .AND. idz.NE.1) then
            DEM(I)%NeighborID(25) = I + LatNx - LatNx*LatNy
        end if
        !  Right Upper Under lattice
        if (idx.NE.LatNx .AND. idy.NE.LatNy .AND. idz.NE.1) then
            DEM(I)%NeighborID(26) = I + 1 + LatNx - LatNx*LatNy
        end if
    end do

    !  Initialize Neighbor Nodelink.
    allocate(IDInner(LatNum))
    allocate(tailInner(LatNum))
    do I = 1,LatNum
        IDInner(I)%No = 0                 !  Array of Lattice's inner particles
        nullify(IDInner(I)%next)
    end do 
#endif

    !  Initialize Periodic DEM Lattice
    if (isPeriodic) then
        write(*,*) "< Parallel Periodic Lattice initializing..."
        LatNum = LatNx*LatNy*LatNz 
        allocate(periodicDEM(LatNum))
        
        do I = 1,LatNum
            !  First element
            periodicDEM(I)%xFlag = 0.0D0
            periodicDEM(I)%yFlag = 0.0D0
            periodicDEM(I)%ID = 0
            nullify(periodicDEM(I)%next)
            !  Tail of the linklist for periodic neighbors
            allocate(Tail)
            Tail%next => periodicDEM(I)
            
            idz = int((I-1)/(LatNx*LatNy))+1
            idy = int((I-(idz-1)*(LatNx*LatNy)-1)/LatNx)+1
            idx = I-(idz-1)*(LatNx*LatNy)-(idy-1)*LatNx
            
            !  Left lattice
            if (idx .EQ. 1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,idy,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Left Upper lattice
            if (idx.EQ.1 .AND. idy.NE.LatNy) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,idy+1,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.1 .AND. idy.EQ.LatNy) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,1,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp                
            else if (idx.NE.1 .AND. idy.EQ.LatNy) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx-1,1,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp                
            end if
            !  Upper lattice
            if (idy .EQ. LatNy) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx,1,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Right Upper lattice
            if (idx.NE.LatNx .AND. idy.EQ.LatNy) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx+1,1,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.LatNx .AND. idy.EQ.LatNy) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,1,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.LatNx .AND. idy.NE.LatNy) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,idy+1,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Left Covered lattice
            if (idx.EQ.1 .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,idy,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Left Upper Covered lattice
            if (idx.EQ.1 .AND. idy.NE.LatNy .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,idy+1,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.1 .AND. idy.EQ.LatNy .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,1,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.NE.1 .AND. idy.EQ.LatNy .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx-1,1,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Upper Covered lattice
            if (idy.EQ.LatNy .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx,1,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Right Upper Covered lattice
            if (idx.NE.LatNx .AND. idy.EQ.LatNy .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx+1,1,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.LatNx .AND. idy.EQ.LatNy .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,1,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.LatNx .AND. idy.NE.LatNy .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,idy+1,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Right Covered lattice
            if (idx.EQ.LatNx .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,idy,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Right Below Covered lattice
            if (idx.EQ.LatNx .AND. idy.NE.1 .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,idy-1,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.LatNx .AND. idy.EQ.1 .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,LatNy,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.NE.LatNx .AND. idy.EQ.1 .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx+1,LatNy,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp               
            end if
            !  Below Covered lattice
            if (idy.EQ.1 .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx,LatNy,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Left Below Covered lattice
            if (idx.EQ.1 .AND. idy.NE.1 .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,idy-1,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.1 .AND. idy.EQ.1 .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,LatNy,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.NE.1 .AND. idy.EQ.1 .AND. idz.NE.LatNz) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx-1,LatNy,idz+1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Right lattice
            if (idx .EQ. LatNx) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,idy,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Right Below lattice
            if (idx.EQ.LatNx .AND. idy.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,idy-1,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.LatNx .AND. idy.EQ.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,LatNy,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.NE.LatNx .AND. idy.EQ.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx+1,LatNy,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Below lattice
            if (idy .EQ. 1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx,LatNy,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Left Below lattice
            if (idx.EQ.1 .AND. idy.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,idy-1,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.1 .AND. idy.EQ.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,LatNy,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.NE.1 .AND. idy.EQ.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx-1,LatNy,idz,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Right Under lattice
            if (idx.EQ.LatNx .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,idy,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Right Below Under lattice
            if (idx.EQ.LatNx .AND. idy.NE.1 .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,idy-1,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.LatNx .AND. idy.EQ.1 .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,LatNy,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.NE.LatNx .AND. idy.EQ.1 .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx+1,LatNy,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Below Under lattice
            if (idy.EQ.1 .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx,LatNy,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Left Below Under lattice
            if (idx.EQ.1 .AND. idy.NE.1 .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,idy-1,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.1 .AND. idy.EQ.1 .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,LatNy,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.NE.1 .AND. idy.EQ.1 .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx-1,LatNy,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,-1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Left Under lattice
            if (idx.EQ.1 .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,idy,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Left Upper Under lattice
            if (idx.EQ.1 .AND. idy.NE.LatNy .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,idy+1,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.1 .AND. idy.EQ.LatNy .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(LatNx,1,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(-1.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.NE.1 .AND. idy.EQ.LatNy .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx-1,1,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Upper Under lattice
            if (idy.EQ.LatNy .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx,1,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
            !  Right Upper Under lattice
            if (idx.EQ.LatNx .AND. idy.NE.LatNy .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,idy+1,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,0.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.EQ.LatNx .AND. idy.EQ.LatNy .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(1,1,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(1.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            else if (idx.NE.LatNx .AND. idy.EQ.LatNy .AND. idz.NE.1) then
                periodicDEM(I)%ID = periodicDEM(I)%ID + 1
                call latticeXY2ID(idx+1,1,idz-1,idLattice,LatNx,LatNy)
                allocate(Temp)
                Temp = PeriodicLattice(0.0D0,1.0D0,idLattice,NULL())
                tail%next%next => Temp
                tail%next => Temp
            end if
        end do
    end if
    
    !!********************** Test of Parallel Lattice ***************************!
    !allocate(test_lattice(3,LatNx*LatNy*LatNz))
    !allocate(test_color(LatNx*LatNy*LatNz))
    !count = 1
    !do K = 1,LatNx
    !    do J = 1,LatNy
    !        do I = 1,LatNx
    !            test_lattice(1,count) = I
    !            test_lattice(2,count) = J
    !            test_lattice(3,count) = K
    !            test_color(count) = 0
    !            count = count + 1
    !        end do
    !    end do
    !end do
    !test_index = 42
    !test_color(test_index) = 2
    !Temp => periodicDEM(test_index)
    !write(*,*) Temp%ID
    !do while(associated(Temp%next))
    !    Temp => Temp%next
    !    write(*,*) Temp%ID
    !    test_color(Temp%ID) = 1
    !end do
    !open(12306,FILE='test_parallel_lattice_2.csv')
    !write(12306,*) 'X,Y,Z,Tag'
    !do I = 1,LatNx*LatNy*LatNz
    !    write(12306,'(F5.1,A2,F5.1,A2,F5.1,A2,I4)') test_lattice(1,I),',',test_lattice(2,I),',',test_lattice(3,I),',',test_color(I)
    !end do
    !stop
    !!********************** Test of Parallel Lattice ***************************!
    
     
    !  Initialize Hertz list. Linklist-Compressed.
    write(*,*) "< Hertz list initializing..."
    
    allocate(Head(NMAX))
    do I = 1,NMAX
        Head(I)%No = 0                    !  Length of Linklist
        !Head(I)%recordTime = 0.0D0        !  record time
        Head(I)%Hertz(1) = 0.0D0          !  Hertz(1)
        Head(I)%Hertz(2) = 0.0D0          !  Hertz(2)
        Head(I)%Hertz(3) = 0.0D0          !  Hertz(3)
        Head(I)%Mrot(1) = 0.0D0           !  Mrot(1)
        Head(I)%Mrot(2) = 0.0D0           !  Mrot(2)
        Head(I)%Mrot(3) = 0.0D0           !  Mrot(3)
        Head(I)%Mtwist(1) = 0.0D0         !  Mtwist(1)
        Head(I)%Mtwist(2) = 0.0D0         !  Mtwist(2)
        Head(I)%Mtwist(3) = 0.0D0         !  Mtwist(3)
        Head(I)%is_touching = .false.     !  Whether touch or not
        Head(I)%is_slipping = .false.     !  Whether slip or not
        Head(I)%is_rolling = .false.      !  Whether roll or not
        Head(I)%is_twisting = .false.     !  Whether twist or not
        nullify(Head(I)%prev)             !  Point to Prev 
        nullify(Head(I)%next)             !  Point to Next
    end do

    return
    end