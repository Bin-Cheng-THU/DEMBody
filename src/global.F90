    !********************************************************************
    !     DEMBody 4.3
    !     ***********
    !
    !     Global parameters.
    !     ------------------
    !
    !********************************************************************

    !#define Grav self_gravity
    
    !#define dims 100
    !#define GravNx 7
    !#define GravNy 7
    !#define GravNz 4

    module global
    implicit none

    !  Parameters
    real(8),parameter :: GravConst = 6.674D-11
    real(8),parameter :: PI = 3.141592653589793D0

    !  Define control parameters of Program
    logical :: isPlanet
    logical :: isRotSystem
    logical :: isQuaternion
    logical :: isContactWall
    logical :: isMovingWall
    logical :: isBondedWall
    logical :: isTriMeshWall
    logical :: isFunnelWall
    logical :: isPeriodic
    logical :: isGravBody
    logical :: isGravTriMesh
    real(8) :: Max_ACC

    !  Conduct linklist
    type :: Nodelink 
        integer :: No
        real(8) :: recordTime
        real(8) :: Hertz(3)
        real(8) :: Mrot(3)
        real(8) :: Mtwist(3)
        logical :: is_touching
        logical :: is_slipping
        logical :: is_rolling
        logical :: is_twisting
        type(Nodelink),pointer :: prev
        type(Nodelink),pointer :: next
    end type Nodelink

    !  Nodelink of Particle
    type(Nodelink),pointer :: Head(:)  
    
    !  Conduct Lattice
    type :: Lattice
        !integer :: ID(3)
#ifdef nConfined
        real(8) :: PositionD(3)
        real(8) :: PositionU(3)
#endif
        integer :: NeighborID(26) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#ifdef self_gravity
!        integer :: GravID
!        real(8) :: Mass
!        real(8) :: MassCenter(3)
!!#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end type Lattice
    
    !  DEM Lattice
    type(Lattice),pointer :: DEM(:)
    
    !  Conduct Neighbor
    type :: Neighbor
        integer :: No
        type(Neighbor),pointer :: next
    end type Neighbor

    !  Nodelink of Neighbor
    type(Neighbor),pointer :: IDInner(:)
    type(Neighbor),pointer :: tailInner(:)

    !  Parallel Lattice params
    real(8) :: LatDx,LatDy,LatDz
    integer :: LatNx,LatNy,LatNz
    integer :: LatNum
    integer,allocatable :: ParallelLatticeColor(:,:)
    
    !  Conduct TriMeshLattice
    type :: GravTriMeshLattice
        real(8) :: PositionCenter(3)
        integer :: NodeID(8)
    end type GravTriMeshLattice
    type(GravTriMeshLattice),pointer :: GravTriMeshGrid(:)
    
    !  Grav Tri Mesh
    real(8),allocatable :: GravTriMeshNode(:,:)

    !  GravTriMeshNode params
    real(8) :: TriLatDx,TriLatDy,TriLatDz
    real(8) :: TriLatMx,TriLatMy,TriLatMz
    integer :: TriLatNx,TriLatNy,TriLatNz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#ifdef self_gravity
!    !  Gravity Lattice params
!    real(8) :: GravDx,GravDy,GravDz
!    integer :: GravNum
!    integer,parameter :: dims = 196
!    integer,parameter :: GravNx = 10
!    integer,parameter :: GravNy = 10
!    integer,parameter :: GravNz = 8
!    
!    !  Conduct Gravity Lattice
!    type :: GravityLattice
!        integer :: num
!        integer :: ID(dims)
!        real(8) :: Mass
!        real(8) :: MassCenter(3)
!    end type
!    
!    type(GravityLattice),pointer :: Gravity(:)
!!#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Define parameters of Particle
    real(8) :: X(3,NMAX),Xdot(3,NMAX),W(3,NMAX)
    real(8) :: Body(NMAX),R(NMAX),Inertia(NMAX)
    real(8) :: F(3,NMAX),FM(3,NMAX)
    real(8) :: Energy(NMAX)
    real(8) :: X0(3,NMAX)
    real(8) :: XT(3,NMAX)
    
    !  Define quaternion of Particle
    real(8) :: Quaternion(4,NMAX)
    
    !  Define material parameters
    real(8) :: m_E,m_nu,m_mu_d,m_mu_s,m_COR
    real(8) :: m_Beta
    real(8) :: m_c,m_r_cut
    real(8) :: m_A
    real(8) :: m_mu_r,m_nita_r
    
    !  Define parameters of Contactable Walls
    integer :: contactWallNum
    integer,allocatable :: contactWallTag(:)
    real(8),allocatable :: contactWallPoint(:,:)
    real(8),allocatable :: contactWallVector(:,:)

    !  Define parameters of Moving Walls
    integer :: movingWallNum
    integer,allocatable :: movingWallTag(:)
    real(8),allocatable :: movingWallPoint(:,:)
    real(8),allocatable :: movingWallVector(:,:)
    
    !  Define parameters of Bonded Walls
    real(8) :: bondedWallX(3),bondedWallXdot(3),bondedWallW(3)
    real(8) :: bondedWallQ(4),bondedWallMatI(3,3),bondedWallMatB(3,3)
    real(8) :: bondedWallBody,bondedWallInertia(3)
    real(8) :: bondedWallF(3),bondedWallFM(3)
    !  body-fixed frame
    real(8) :: bondedWallWB(3),bondedWallWdotB(3)    
    real(8),allocatable :: bondedWallMeshPoint(:,:) 
    integer :: bondedWallNum
    integer,allocatable :: bondedWallTag(:)
    real(8),allocatable :: bondedWallPoint(:,:)
    real(8),allocatable :: bondedWallVectorN(:,:)
    real(8),allocatable :: bondedWallVectorTx(:,:)
    real(8),allocatable :: bondedWallVectorTy(:,:)
    real(8),allocatable :: bondedWallLx(:),bondedWallLy(:)
    
    !  Define parameters of TriMesh Walls
    integer :: trimeshWallNum
    integer,allocatable :: trimeshWallTag(:)
    real(8),allocatable :: trimeshWallPoint(:,:)
    real(8),allocatable :: trimeshWallVectorN(:,:)
    real(8),allocatable :: trimeshWallVectorTx(:,:)
    real(8),allocatable :: trimeshWallVectorTy(:,:)
    real(8),allocatable :: trimeshWallLength(:,:)  !  Tx2, Tx*Ty, Ty2, L
    
    !  Define parameters of Funnel Walls
    integer :: funnelWallNum
    integer,allocatable :: funnelWallTag(:)
    real(8),allocatable :: funnelWallPoint(:,:)
    real(8),allocatable :: funnelWallVector(:,:)
    real(8),allocatable :: funnelWallRadius(:,:)
    real(8),allocatable :: funnelWallLength(:)
    
    !  Define parameters of Periodic Walls
    real(8) :: PlaSx1p(3),PlaSx1v(3)   !  Position of WallX1
    real(8) :: PlaSx2p(3),PlaSx2v(3)   !  Position of WallX2
    real(8) :: PlaSy1p(3),PlaSy1v(3)   !  Position of WallY1
    real(8) :: PlaSy2p(3),PlaSy2v(3)   !  Position of WallY2
    real(8) :: LenBoxX
    real(8) :: LenBoxY
    real(8) :: gamma
    integer :: Tag1(NMAX),Tag2(NMAX),Tag3(NMAX),Tag4(NMAX)
    
    !Define parameters of GravBody
    integer :: gravBodyTag
    real(8) :: gravBodyX(3),gravBodyXdot(3),gravBodyW(3)
    real(8) :: gravBodyQ(4)
    real(8) :: gravBodyBody,gravBodyR,gravBodyInertia
    real(8) :: gravBodyF(3),gravBodyFM(3)
        
    !  Define parameters of Saturn and Pan
    real(8) :: muS
    real(8) :: muP
    real(8) :: rOrig(3)
    real(8) :: omega
    
    !  Define parameters of Rotary system
    real(8) :: sysOmega
    real(8) :: sysGrav
    
    !  Define parameters of Gravity
    real(8) :: G
    
    !  Define parameters of Mesh
    integer :: N
    real(8) :: Dx,Dy,Dz
    real(8) :: Mx,My,Mz
    integer :: Nx,Ny,Nz
    integer :: Linklist(NMAX)
    real(8) :: verlet

    !  Define parameters of Program
    real(8) :: T1,T2
    real(8) :: Time,Deltat,Tcrit,Dt,Tnext
    real(8) :: CheckPointDt, CheckPointTnext
    integer :: Step
    logical :: refreshLattice
    integer :: refreshNum
   
    end module