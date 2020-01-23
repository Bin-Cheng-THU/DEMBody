    !********************************************************************
    !     DEMBody 8.3
    !     ***********
    !
    !     Global parameters.
    !     ------------------
    !
    !********************************************************************

    module global
    implicit none

    !DEM: 存储网格的相邻网格编号，采用数组形式存储，共LatNum个
    !periodicDEM: 存储网格的“周期相邻”网格编号，采用链表形式存储，共LatNum个
    !!!!!mixDEM：存储网格的相邻网格编号，采用链表形式存储所有网格（包括自身网格），本网格为适用于大颗粒的粗网格

    !IDInner： 存储网格内的颗粒编号，共LatNum个
    !trimeshDEM： 存储网格内的三角面元编号，共LatNum个
    !bondedTriMeshDEM： 存储网格内的bonded三角面元编号，共LatNum个
    !biDisperseDEM： 存储网格内的biDisperse颗粒编号，共LatNum个
    !!!!!MixIDInner 存储网格内的BiDisperse颗粒编号，注意此网格为用于biParticle的粗网格
    ![initialXXX: initial the DEM lattice; only in the initial step;
    !latticeGenerateXXX: generate the lattice structure, which stores the ID of the bodies; every time step]

    !Head： 存储颗粒接触编号及接触历史信息，共NMAX个
    !HeadBiDisperse： 存储biDisperse颗粒接触编号及接触历史信息，共biDisperseNum个[颗粒与biDisperse的接触存储与Head中！！！]
    ![注意虽然检索颗粒时存在团聚搜索的结构，即每个网格内保序，各个网格按照顺序，但网格之间的颗粒序号并不是保序的；但在存储接触历史时均按照保序操作，这样使得在颗粒穿过网格时查找接触历史的操作仍能够正常进行。同时注意，镜像颗粒与颗粒的接触存储是统一的，即均按照颗粒序号排序。]

    !GravTriMeshGrid： 存储三角元对应网格的质心与网格点编号


    !  Parameters
    real(8),parameter :: GravConst = 6.674184D-11
    real(8),parameter :: PI = 3.141592653589793D0
    character(10),parameter :: VERSION = '8.3'

    !  Define control parameters of Program
    character(10) :: vsDEMBody
    logical :: isPlanet
    logical :: isRotSystem
    logical :: isQuaternion
    logical :: isContactWall
    logical :: isMovingWall
    logical :: isBondedWall
    logical :: isTriMeshWall
    logical :: isBondedTriMeshWall
    logical :: isFunnelWall
    logical :: isPeriodic
    logical :: isGravBody
    logical :: isSphereBody
    logical :: isBiDisperse
    logical :: isGravTriMesh
    real(8) :: Max_ACC

    !  Define parameters of Program
    real(8) :: T1,T2
    real(8) :: Time,Deltat,Tcrit,Dt,Tnext
    real(8) :: CheckPointDt, CheckPointTnext
    integer :: Step
    logical :: refreshLattice
    integer :: refreshNum

    !  Conduct linklist
    type :: Nodelink 
        integer :: No
        !real(8) :: recordTime
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
    
    !  Define parameters of Parallel Lattice
    integer :: N
    real(8) :: LatDx,LatDy,LatDz
    integer :: LatNx,LatNy,LatNz
    real(8) :: LatMx,LatMy,LatMz
    integer :: LatNum
    integer :: Linklist(NMAX)
    real(8) :: verlet
    !integer,allocatable :: ParallelLatticeColor(:,:)
    
    !  Conduct Lattice
    type :: Lattice
        !integer :: ID(3)
#ifdef nConfined
        real(8) :: PositionD(3)
        real(8) :: PositionU(3)
#endif
        integer :: NeighborID(26) 
    end type Lattice
    !  DEM Lattice
    type(Lattice),pointer :: DEM(:)
    
    !  Periodic Lattice
    type :: PeriodicLattice
        real(8) :: xFlag
        real(8) :: yFlag
        integer :: ID
        type(PeriodicLattice),pointer :: next
    end type PeriodicLattice
    !  DEM Periodic Lattice
    type(PeriodicLattice),pointer :: periodicDEM(:)
        
    !  Conduct Neighbor
    type :: Neighbor
        integer :: No
        type(Neighbor),pointer :: next
    end type Neighbor
    !  Nodelink of Neighbor
    type(Neighbor),pointer :: IDInner(:)
    type(Neighbor),pointer :: tailInner(:)  
    
    !  GravTriMeshNode params
    real(8) :: TriLatDx,TriLatDy,TriLatDz
    real(8) :: TriLatMx,TriLatMy,TriLatMz
    integer :: TriLatNx,TriLatNy,TriLatNz    
    !  Conduct TriMeshLattice
    type :: GravTriMeshLattice
        real(8) :: PositionCenter(3)
        integer :: NodeID(8)
    end type GravTriMeshLattice
    type(GravTriMeshLattice),pointer :: GravTriMeshGrid(:)    
    !  Grav Tri Mesh
    real(8),allocatable :: GravTriMeshNode(:,:)

    !  Define parameters of Particle
    real(8) :: X(3,NMAX),Xdot(3,NMAX),W(3,NMAX)
    real(8) :: Body(NMAX),R(NMAX),Inertia(NMAX)
    real(8) :: F(3,NMAX),FM(3,NMAX)
    real(8) :: Energy(NMAX)
    real(8) :: Heat(3,NMAX)
    real(8) :: X0(3,NMAX)
    real(8) :: XT(3,NMAX)
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
    real(8),allocatable :: movingWallPointInit(:,:)
    real(8),allocatable :: movingWallVector(:,:)    
    real(8),allocatable :: movingWallVelocity(:,:)
    real(8) :: movingWallTstart,movingWallTend
    real(8) :: movingWallA,movingWallOmega
    real(8) :: movingWallNormal(3)
    integer :: currentStep
    !real(8),allocatable :: movingWallPointStore(:,:,:)
    !real(8),allocatable :: movingWallVectorStore(:,:,:)
    !real(8),allocatable :: movingWallVelocityStore(:,:,:)
    
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
    !  Define trimesh in lattice
    type :: trimeshLattice
        integer :: No
        type(trimeshLattice),pointer :: next
    end type trimeshLattice
    type(trimeshLattice),pointer :: trimeshDEM(:)
    
    !  Define parameters of Bonded TriMesh Walls
    real(8) :: bondedTriMeshWallX(3),bondedTriMeshWallXdot(3),bondedTriMeshWallW(3)
    real(8) :: bondedTriMeshWallQ(4),bondedTriMeshWallMatI(3,3),bondedTriMeshWallMatB(3,3)
    real(8) :: bondedTriMeshWallBody,bondedTriMeshWallInertia(3)
    real(8) :: bondedTriMeshWallF(3),bondedTriMeshWallFM(3)
    !  body-fixed frame
    integer :: bondedTriMeshWallNum
    integer,allocatable :: bondedTriMeshWallTag(:)
    real(8),allocatable :: bondedTriMeshWallPoint(:,:)
    real(8),allocatable :: bondedTriMeshWallVectorN(:,:)
    real(8),allocatable :: bondedTriMeshWallVectorTx(:,:)
    real(8),allocatable :: bondedTriMeshWallVectorTy(:,:)
    real(8),allocatable :: bondedTriMeshWallLength(:,:)  !  Tx2, Tx*Ty, Ty2, L
    type(triMeshLattice),pointer :: bondedTriMeshDEM(:)
    
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
    
    !  Define parameters of GravBody
    integer :: gravBodyTag
    real(8) :: gravBodyX(3),gravBodyXdot(3),gravBodyW(3)
    real(8) :: gravBodyQ(4)
    real(8) :: gravBodyBody,gravBodyR,gravBodyInertia
    real(8) :: gravBodyF(3),gravBodyFM(3)
    
    !  Define parameters of SphereBody
    integer :: sphereBodyNum
    integer,allocatable :: sphereBodyTag(:)
    real(8),allocatable :: sphereBodyX(:,:),sphereBodyXdot(:,:),sphereBodyW(:,:)
    real(8),allocatable :: sphereBodyQ(:,:)
    real(8),allocatable :: sphereBodyBody(:),sphereBodyR(:),sphereBodyInertia(:)
    real(8),allocatable :: sphereBodyF(:,:),sphereBodyFM(:,:)

    !******************************************** v8.0 ***********************************************
    !  Define material properties for BiDisperse [maybe we could add other properties]
    real(8) :: m_Ei

    !  Define parameters of BiDisperse
    integer :: biDisperseNum
    integer,allocatable :: biDisperseTag(:)  !  Unique tag for biDisperse particles, numbered with walls and particles. For particle/BiDisperse particle contact
    real(8),allocatable :: biDisperseX(:,:),biDisperseXdot(:,:),biDisperseW(:,:)
    real(8),allocatable :: biDisperseQ(:,:)
    real(8),allocatable :: biDisperseBody(:),biDisperseR(:),biDisperseInertia(:)
    real(8),allocatable :: biDisperseF(:,:),biDisperseFM(:,:)
    real(8),allocatable :: biDisperseXT(:,:)
    real(8),allocatable :: biDisperseEnergy(:)
    real(8),allocatable :: biDisperseHeat(:)  !  Only slip energy

    !  Nodelink of BiDisperse
    type(Nodelink),pointer :: HeadBiDisperse(:)  !  For contact between BiDisperse particles

    !  Define structure for particle-biParticle contact
    !integer :: biDisperseScale  !  biDisperseR/LatDx == biDisperseR/Rmax/2.5
    !  Define bidisperse in lattice
    type :: biDisperseLattice
        integer :: No
        type(biDisperseLattice),pointer :: next
    end type biDisperseLattice
    type(biDisperseLattice),pointer :: biDisperseDEM(:)
    type(biDisperseLattice),pointer :: biDisperseDEMtail(:)

    !  Define refresh structure, inheriting from MixBiMesh
    logical :: refreshBiDisperseLattice
    integer :: refreshBiDisperseNum
    real(8) :: verletBiDisperse  !  (LatDx-Rmax)/2 = (2.5Rmax-Rmax)/2
    
    !******************************************** v8.0 ***********************************************
        
    !  Define parameters of Saturn and Pan
    real(8) :: muS
    real(8) :: muP
    real(8) :: rOrig(3)
    real(8) :: omega
    
    !  Define parameters of Rotary system
    real(8) :: sysOmega

    !  Define parameters of YORP evolution in local Frame
    logical :: islocalYORP
    real(8) :: localYORPTstart,localYORPTend
    real(8) :: localYORPOmega
    real(8) :: localYORPIncrement,localYORPDt
    real(8) :: localYORPTnext
    real(8) :: localYORPdensity
    !******************************************** v8.3 ***********************************************
    !  biDisperse Particles + local YORP
    real(8) :: localYORPcenter(3)
    real(8) :: localYORPradius
    real(8) :: localYORPmass
    
    !  Define parameters of Gravity
    real(8) :: G(3)
   
    end module