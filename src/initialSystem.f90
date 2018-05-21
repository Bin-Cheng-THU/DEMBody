    !********************************************************************
    !     DEMBody 4.0
    !     ***********
    !
    !     Initialization of global scalars.
    !     ---------------------------------
    !
    !********************************************************************  
    SUBROUTINE initialSystem()

    use global
    use loadFile
    implicit none
    
    integer I,J,K
    integer :: wallFlag
    integer :: nRow
    integer :: idx,idy,idz
    integer :: idgx,idgy,idgz
    real(8) :: lx,ly,lz

    !  Initialize parameters and set useful constants.
    Time = 0.0D0
    Tnext = 0.0D0
    Step = 0
    CheckPointTnext = 0.0D0
    
    !  Initialize control parameters
    write(*,*) "System file loading..."
    
    open (1000,FILE="../Input/systemControl.txt",STATUS='OLD',BLANK='NULL',POSITION='REWIND')

    read (1000,*)                      !  DEM version
    read (1000,*) 
    read (1000,*) LatDx,LatDy,LatDz    !  Parallel Lattice grid interval
    read (1000,*) LatNx,LatNy,LatNz    !  Parallel Lattice grid number
    read (1000,*)
    read (1000,*) Deltat,Tcrit         !  output time interval and termination time
    read (1000,*) CheckPointDt         !  check point time interval
    read (1000,*) Dt                   !  time step of second-step integral
    read (1000,*)                         
    read (1000,*) Dx,Dy,Dz             !  mesh grid interval
    read (1000,*) Mx,My,Mz             !  mesh grid origin
    read (1000,*) Nx,Ny,Nz             !  mesh grid number
    read (1000,*)                         
    read (1000,*) isPlanet             !  whether use Planet Gravity function
    read (1000,*) isQuaternion         !  whether intergrate Quaternion
    read (1000,*) isContactWall        !  whether use Contactable Walls
    read (1000,*) isMovingWall         !  whether use Moving Walls
    read (1000,*) isBondedWall         !  whether use Bonded Walls
    read (1000,*) isFunnelWall         !  whether use Funnel Walls
    read (1000,*) isPeriodic           !  whether use Periodic function
    read (1000,*) isGravBody           !  whether use Gravity Body
    read (1000,*) MAX_ACC              !  maximum contact acceleration
    read (1000,*)                         
    read (1000,*) G                    !  the totle gravity    
    
    wallFlag = 1
    
    !  initial the contactable walls
    if (isContactWall) then
        write(*,*) 'is Contactable walls, loading...'
        read (1000,*) 
        read (1000,*) contactWallNum
        allocate (contactWallTag(contactWallNum))
        allocate (contactWallPoint(3,contactWallNum))
        allocate (contactWallVector(3,contactWallNum))
        do I = 1,contactWallNum
            contactWallTag(I) = N + wallFlag
            wallFlag = wallFlag + 1
            read (1000,*) (contactWallPoint(K,I),K=1,3),(contactWallVector(K,I),K=1,3)
        end do
    else
        read (1000,*) 
        read (1000,*)
    end if
    
    !  initial the bonded walls
    if (isBondedWall) then
        write(*,*) 'is Bonded walls, loading...'
        read (1000,*)
        read (1000,*) (bondedWallX(K),K=1,3),(bondedWallXdot(K),K=1,3),(bondedWallW(K),K=1,3)
        read (1000,*) (bondedWallQ(K),K=1,4)
        read (1000,*) bondedWallBody,(bondedWallInertia(K),K=1,3)
        read (1000,*) bondedWallNum
        allocate (bondedWallTag(bondedWallNum))
        allocate (bondedWallPoint(3,bondedWallNum))
        allocate (bondedWallVectorN(3,bondedWallNum))
        allocate (bondedWallVectorTx(3,bondedWallNum))
        allocate (bondedWallVectorTy(3,bondedWallNum))
        allocate (bondedWallLx(bondedWallNum))
        allocate (bondedWallLy(bondedWallNum))
        do I = 1,bondedWallNum
            bondedWallTag(I) = N + wallFlag
            wallFlag = wallFlag + 1
            read (1000,*) (bondedWallPoint(K,I),K=1,3),(bondedWallVectorN(K,I),K=1,3),(bondedWallVectorTx(K,I),K=1,3),(bondedWallVectorTy(K,I),K=1,3),bondedWallLx(I),bondedWallLy(I)
        end do
        call attitudeBondedWallsQM
        !  unified input parameters
        do K = 1,3
            bondedWallWB(K) = bondedWallMatI(K,1)*bondedWallW(1) + bondedWallMatI(K,2)*bondedWallW(2) + bondedWallMatI(K,3)*bondedWallW(3)
        end do
        !  load bonded wall mesh
        open (1500,File="../Input/bondedWallPoint.vtk")
        nRow = GetFileN(1500)
        allocate (bondedWallMeshPoint(3,nRow))
        do I = 1,nRow
            read (1500,*) (bondedWallMeshPoint(K,I),K=1,3)
        end do
        close(1500)
    else
        read (1000,*)
        read (1000,*)
    end if
    
    !  initial the funnel walls
    if (isFunnelWall) then
        write(*,*) 'is Funnel walls, loading...'
        read (1000,*) 
        read (1000,*) funnelWallNum
        allocate (funnelWallTag(funnelWallNum))
        allocate (funnelWallPoint(3,funnelWallNum))
        allocate (funnelWallVector(3,funnelWallNum))
        allocate (funnelWallRadius(2,funnelWallNum))
        allocate (funnelWallLength(funnelWallNum))
        do I = 1,funnelWallNum
            funnelWallTag(I) = N + wallFlag
            wallFlag = wallFlag + 1
            read (1000,*) (funnelWallPoint(K,I),K=1,3),(funnelWallVector(K,I),K=1,3),(funnelWallRadius(K,I),K=1,2),funnelWallLength(I)
        end do
    else
        read (1000,*)
        read (1000,*)
    end if
    
    !  initial the geometry boundary of the Particle Box
    if (isPeriodic) then
        write(*,*) 'is Periodic boundary, loading...'
        read (1000,*)
        read (1000,*) (PlaSx1p(K),K=1,3),(PlaSx1v(K),K=1,3)
        read (1000,*) (PlaSx2p(K),K=1,3),(PlaSx2v(K),K=1,3)
        read (1000,*) (PlaSy1p(K),K=1,3),(PlaSy1v(K),K=1,3)
        read (1000,*) (PlaSy2p(K),K=1,3),(PlaSy2v(K),K=1,3)
        read (1000,*) LenBoxX
        read (1000,*) LenBoxY
        read (1000,*) gamma
    else
        read (1000,*)
        read (1000,*)
    end if
    
    !  initial the gravity body
    if (isGravBody) then
        write(*,*) 'is Gravity Body, loading...'
        read (1000,*)
        gravBodyTag = N + wallFlag
        read (1000,*) (gravBodyX(K),K=1,3),(gravBodyXdot(K),K=1,3),(gravBodyW(K),K=1,3)
        read (1000,*) gravBodyBody,gravBodyR,gravBodyInertia
        gravBodyQ(1) = 0.0D0
        gravBodyQ(2) = 0.0D0
        gravBodyQ(3) = 0.0D0
        gravBodyQ(4) = 1.0D0
    else
        read (1000,*)
        read (1000,*)
    end if
    close(1000)
    
    !  read properties of Saturn and Pan
    if (isPlanet) then
        write(*,*) 'is Planet system, loading...'
        open (2000,FILE="../Input/planetProperties.txt",STATUS='OLD',BLANK='NULL',POSITION='REWIND')
        read (2000,*)
        read (2000,*)  muS                 !  gravity constant of Saturn
        read (2000,*)
        read (2000,*)  muP                 !  gravity constant of Pan
        read (2000,*)
        read (2000,*)  omiga               !  angular velocity of Pan's initial revolution (suppose tidal locking and anticlockwise, usually negative)
        read (2000,*)
        read (2000,*)  (rOrig(K),K=1,3)    !  Saturn to Pan's initial position, i.e., the origin point
        close(2000)
    end if
    
!!#ifdef self_gravity
!    !  Initialize Gravity Lattice.
!    write(*,*) "Gravity Lattice initializing..."
!    
!    GravNum = GravNx*GravNy*GravNz
!    GravDx = LatDx*(LatNx/GravNx)
!    GravDy = LatDy*(LatNy/GravNy)
!    GravDz = LatDz*(LatNz/GravNz)
!
!    allocate(Gravity(GravNum))
!    allocate(ParallelLatticeColor(1,NMAX))
!    do I = 1,GravNum
!        Gravity(I)%num = 0
!        Gravity(I)%ID = 0
!        Gravity(I)%Mass = 0.0D0
!        Gravity(I)%MassCenter = 0.0D0
!    end do
!!#end if

    !  Initialize Parallel Lattice.
    write(*,*) "Parallel Lattice initializing..."
    
    LatNum = LatNx*LatNy*LatNz             !  Paraller Number
    allocate(DEM(LatNum))
    do I = 1,LatNum
        idz = int((I-1)/(LatNx*LatNy))+1
        idy = int((I-(idz-1)*(LatNx*LatNy)-1)/LatNx)+1
        idx = I-(idz-1)*(LatNx*LatNy)-(idy-1)*LatNx
        lx = dble(idx-1)*LatDx - Mx
        ly = dble(idy-1)*LatDy - My
        lz = dble(idz-1)*LatDz - Mz
        DEM(I)%ID(1) = idx                 !  IDx of Lattice
        DEM(I)%ID(2) = idy                 !  IDy of Lattice
        DEM(I)%ID(3) = idz                 !  IDz of Lattice
        DEM(I)%PositionD(1) = lx           !  PositionDx of Lattice
        DEM(I)%PositionD(2) = ly           !  PositionDy of Lattice
        DEM(I)%PositionD(3) = lz           !  PositionDz of Lattice
        DEM(I)%PositionU(1) = lx + LatDx   !  PositionUx of Lattice
        DEM(I)%PositionU(2) = ly + LatDy   !  PositionUy of Lattice
        DEM(I)%PositionU(3) = lz + LatDz   !  PositionUz of Lattice        
        DEM(I)%NoInner = 0                 !  Length of Lattice's inner particles
        DEM(I)%NoOuter = 0                 !  Length of Lattice's outer particles
        DEM(I)%IDInner = 0                 !  Array of Lattice's inner particles
        DEM(I)%IDOuter = 0                 !  Array of lattice's outer particles
!!#ifdef self_gravity
!        idgx = int((idx-1)/(LatNx/GravNx))+1
!        idgy = int((idy-1)/(LatNy/GravNy))+1
!        idgz = int((idz-1)/(LatNz/GravNz))+1
!        DEM(I)%GravID = idgx + (idgy-1)*GravNx + (idgz-1)*GravNx*GravNy
!        DEM(I)%Mass = 0.0D0
!        DEM(I)%MassCenter = DEM(I)%PositionD
!        
!        Gravity(DEM(I)%GravID)%num = Gravity(DEM(I)%GravID)%num + 1
!        Gravity(DEM(I)%GravID)%ID(Gravity(DEM(I)%GravID)%num) = I
!!#endif
    end do

    !  Initialize Hertz list. Linklist-Compressed.
    write(*,*) "Hertz list initializing..."
    
    allocate(Head(NMAX))
    do I = 1,NMAX
        Head(I)%No = 0                    !  Length of Linklist
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