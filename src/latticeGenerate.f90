    !********************************************************************
    !     DEMBody 3.0
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
    !********************************************************************
    subroutine latticeGenerate()

    use global
    implicit none

    integer :: I,J
    integer :: Tag
    logical :: Flag
    real(8) :: positionD(3),positionU(3)
    integer :: color
    
    !  initialized
    do I = 1,LatNum
        
        !  reset Inner & Outer quantities
        DEM(I)%NoInner = 0               
        DEM(I)%NoOuter = 0               
        DEM(I)%IDInner = 0               
        DEM(I)%IDOuter = 0               
    
    end do
    
    do I = 1,N
        
        Tag = floor((X(1,I)+Mx)/LatDx) + 1 + floor((X(2,I)+My)/LatDy)*LatNx + floor((X(3,I)+Mz)/LatDz)*(LatNx*LatNy)
        
        !****************
        !Note that if the particle is out of computational region, the Tag
        !may be beyond Normal Parallel Lattice, which would generat error when put into
        !DEM array; so we must check the position of the particle.
        !****************
        
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
    end do

    !****************
    !Test Parallel Lattice Method using color visualization.
    !****************
            
    end