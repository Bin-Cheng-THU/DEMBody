    !********************************************************************
    !     DEMBody 7.0
    !     ***********
    !
    !     Initialization of biDisperse.
    !     ---------------------------------
    !
    !********************************************************************
    SUBROUTINE initialBiDisperse()

    use global
    implicit none

    integer I,J,K,L
    integer :: idx,idy,idz
    
    integer :: idLattice
    type(MixLattice),pointer :: Temp
    type(MixLattice),pointer :: Tail
    real(8) :: ostart,oend

    refreshBiDisperseLattice = .true.
    refreshBiDisperseNum = 0    

    !  initialize biDisperse
    LatNum = LatNx*LatNy*LatNz 
    allocate(biDisperseDEM(LatNum))
    allocate(biDisperseDEMtail(LatNum))
    do I = 1,LatNum
        biDisperseDEM(I)%No = 0
        nullify(biDisperseDEM(I)%next)
    end do
    
    !  initialize  MixMesh
    MixLatNum = MixLatNx*MixLatNy*MixLatNz
    allocate(MixIDInner(MixLatNum))
    allocate(MixtailInner(MixLatNum))
    do I = 1,MixLatNum
        MixIDInner(I)%No = 0                 !  Array of MixMeshLattice's inner particles
        nullify(MixIDInner(I)%next)
    end do 

    !  initialize  HeadBiDisperse
    allocate(HeadBiDisperse(biDisperseNum))
    do I = 1,biDisperseNum
        HeadBiDisperse(I)%No = 0                    !  Length of Linklist
        !HeadBiDisperse(I)%recordTime = 0.0D0        !  record time
        HeadBiDisperse(I)%Hertz(1) = 0.0D0          !  Hertz(1)
        HeadBiDisperse(I)%Hertz(2) = 0.0D0          !  Hertz(2)
        HeadBiDisperse(I)%Hertz(3) = 0.0D0          !  Hertz(3)
        HeadBiDisperse(I)%Mrot(1) = 0.0D0           !  Mrot(1)
        HeadBiDisperse(I)%Mrot(2) = 0.0D0           !  Mrot(2)
        HeadBiDisperse(I)%Mrot(3) = 0.0D0           !  Mrot(3)
        HeadBiDisperse(I)%Mtwist(1) = 0.0D0         !  Mtwist(1)
        HeadBiDisperse(I)%Mtwist(2) = 0.0D0         !  Mtwist(2)
        HeadBiDisperse(I)%Mtwist(3) = 0.0D0         !  Mtwist(3)
        HeadBiDisperse(I)%is_touching = .false.     !  Whether touch or not
        HeadBiDisperse(I)%is_slipping = .false.     !  Whether slip or not
        HeadBiDisperse(I)%is_rolling = .false.      !  Whether roll or not
        HeadBiDisperse(I)%is_twisting = .false.     !  Whether twist or not
        nullify(HeadBiDisperse(I)%prev)             !  Point to Prev 
        nullify(HeadBiDisperse(I)%next)             !  Point to Next
    end do

    !  initialize mixDEM
    allocate(mixDEM(MixLatNum))
    do I = 1,MixLatNum
        !  First element
        mixDEM(I)%NeighborID = 0
        nullify(mixDEM(I)%next)
        !  Tail of the linklist for periodic neighbors
        allocate(Tail)
        Tail%next => mixDEM(I)

        idz = int((I-1)/(MixLatNx*MixLatNy))+1
        idy = int((I-(idz-1)*(MixLatNx*MixLatNy)-1)/MixLatNx)+1
        idx = I-(idz-1)*(MixLatNx*MixLatNy)-(idy-1)*MixLatNx

         !  Center lattice
        if (.True.) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Left lattice
        if (idx .NE. 1) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I-1
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Left Upper lattice
        if (idx.NE.1 .AND. idy.NE.MixLatNy) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I - 1 + MixLatNx
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Upper lattice
        if (idy .NE. MixLatNy) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I + MixLatNx
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Right Upper lattice
        if (idx.NE.MixLatNx .AND. idy.NE.MixLatNy) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I + 1 + MixLatNx
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Center Covered lattice
        if (idz .NE. MixLatNz) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I + MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Left Covered lattice
        if (idx.NE.1 .AND. idz.NE.MixLatNz) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I -1 + MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Left Upper Covered lattice
        if (idx.NE.1 .AND. idy.NE.MixLatNy .AND. idz.NE.MixLatNz) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I - 1 + MixLatNx + MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Upper Covered lattice
        if (idy.NE.MixLatNy .AND. idz.NE.MixLatNz) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I + MixLatNx + MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Right Upper Covered lattice
        if (idx.NE.MixLatNx .AND. idy.NE.MixLatNy .AND. idz.NE.MixLatNz) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I + 1 + MixLatNx + MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Right Covered lattice
        if (idx.NE.MixLatNx .AND. idz.NE.MixLatNz) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I + 1 + MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Right Below Covered lattice
        if (idx.NE.MixLatNx .AND. idy.NE.1 .AND. idz.NE.MixLatNz) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I + 1 - MixLatNx + MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Below Covered lattice
        if (idy.NE.1 .AND. idz.NE.MixLatNz) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I - MixLatNx + MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Left Below Covered lattice
        if (idx.NE.1 .AND. idy.NE.1 .AND. idz.NE.MixLatNz) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I - 1 - MixLatNx + MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Right lattice
        if (idx .NE. MixLatNx) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I + 1
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Right Below lattice
        if (idx.NE.MixLatNx .AND. idy.NE.1) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I + 1 - MixLatNx
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Below lattice
        if (idy .NE. 1) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I - MixLatNx
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Left Below lattice
        if (idx.NE.1 .AND. idy.NE.1) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I - 1 - MixLatNx
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Center Under lattice
        if (idz .NE. 1) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I - MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Right Under lattice
        if (idx.NE.MixLatNx .AND. idz.NE.1) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I + 1 - MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Right Below Under lattice
        if (idx.NE.MixLatNx .AND. idy.NE.1 .AND. idz.NE.1) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I + 1 - MixLatNx - MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Below Under lattice
        if (idy.NE.1 .AND. idz.NE.1) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I - MixLatNx - MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Left Below Under lattice
        if (idx.NE.1 .AND. idy.NE.1 .AND. idz.NE.1) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I - 1 - MixLatNx - MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Left Under lattice
        if (idx.NE.1 .AND. idz.NE.1) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I - 1 - MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Left Upper Under lattice
        if (idx.NE.1 .AND. idy.NE.MixLatNy .AND. idz.NE.1) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I - 1 + MixLatNx - MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Upper Under lattice
        if (idy.NE.MixLatNy .AND. idz.NE.1) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I + MixLatNx - MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
        !  Right Upper Under lattice
        if (idx.NE.MixLatNx .AND. idy.NE.MixLatNy .AND. idz.NE.1) then
            mixDEM(I)%NeighborID = mixDEM(I)%NeighborID + 1
            idLattice = I + 1 + MixLatNx - MixLatNx*MixLatNy
            allocate(Temp)
            Temp = MixLattice(idLattice,NULL())
            tail%next%next => Temp
            tail%next => Temp
        end if
    end do

    end