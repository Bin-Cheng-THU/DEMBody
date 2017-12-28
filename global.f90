    !********************************************************************
    !     DEMbody 2.0
    !     ***********
    !
    !     Global parameters.
    !     ------------------
    !
    !********************************************************************
    module global
    implicit none

    !  Parameters
    integer,parameter :: NMAX = 300
    real(8),parameter :: PI = 3.141592653589793

    !  Define control parameters of Program
    logical :: isPlanet
    logical :: isQuaternion
    logical :: isContactWall
    logical :: isMovingWall
    logical :: isFunnelWall
    logical :: isPeriodic
    real(8) :: Max_ACC

    !  Conduct linklist
    type :: Nodelink 
        integer :: No
        real    :: Hertz(3)
        real    :: Mrot(3)
        logical :: is_touching
        logical :: is_slipping
        logical :: is_rolling
        type(Nodelink),pointer :: prev
        type(Nodelink),pointer :: next
    end type Nodelink

    !  Nodelink of Particle
    type(Nodelink),pointer :: Head(:)  
    
    !  Define parameters of Particle
    real(8) :: X(3,NMAX),Xdot(3,NMAX),W(3,NMAX)
    real(8) :: Body(NMAX),R(NMAX),Inertia(NMAX)
    real(8) :: F(3,NMAX),FM(3,NMAX)
    real(8) :: Energy(NMAX)
    real(8) :: X0(3,NMAX)
    
    !  Define quaternion of Particle
    real(8) :: Quaternion(4,NMAX)
    
    !  Define material parameters
    real(8) :: m_E,m_nu,m_mu_d,m_mu_s,m_COR,m_mu_r,m_nita_r,m_Beta,m_r_cut,m_A
    
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
    real(8) :: LenBox
    real(8) :: gamma
    integer :: Tag1(NMAX),Tag2(NMAX),Tag3(NMAX),Tag4(NMAX)
    
    !  Define parameters of Saturn and Pan
    real(8) :: muS
    real(8) :: muP
    real(8) :: rOrig(3)
    real(8) :: omiga
    
    !  Define parameters of Gravity
    real(8) :: G
    
    !  Define parameters of Mesh
    integer :: N
    real(8) :: Dx,Dy,Dz
    real(8) :: Mx,My,Mz
    integer :: Nx,Ny,Nz
    integer :: Linklist(NMAX)

    !  Define parameters of Program
    real(8) :: T1,T2
    real(8) :: Time,Deltat,Tcrit,Dt,Tnext
    real(8) :: CheckPointDt, CheckPointTnext
    integer :: Step
   
    end module
    
    !*
    !历史遗留问题
    !COMMON/NBODY/  X,Xdot,W,Body,R,Inertia,F,FM,Energy,Head,X0 ! Particle - Particle
    !!
    !COMMON/MATERIAL/ m_E,m_nu,m_mu_d,m_mu_s,m_COR,m_mu_r,m_nita_r,m_Beta,m_r_cut,m_A,G !Material param
    !!
    !COMMON/PROPERTIES/ muS,muP,rOrig,omiga  !  Dynamical system param
    !!
    !COMMON/WALLS/  PlaSx1p,PlaSx1v,PlaSx2p,PlaSx2v,PlaSy1p,PlaSy1v,PlaSy2p,PlaSy2v,& ! Particle - Wall  
    !&LenBox,gamma,Tag1,Tag2,Tag3,Tag4 !Do not need to add the contact walls because these parameters are in Module
    !!
    !COMMON/MESH/  N,Dx,Dy,Dz,Mx,My,Mz,Nx,Ny,Nz,Linklist  ! Mesh param
    !!
    !COMMON/PARAMS/ T1,T2,Time,Deltat,Tcrit,Dt,Tnext,Step,CheckPointDt,CheckPointTnext ! Console param
    !!
    !COMMON/CONTROL/ isGravity,isPlanet,isPeriodic,Max_ACC  !  Control param   
!* 