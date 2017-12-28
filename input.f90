    !********************************************************************
    !     DEMbody 2.0
    !     ***********
    !
    !     Parameter input.
    !     ----------------
    !
    !********************************************************************
    subroutine input()

    use global
    implicit none
    
    integer :: I,K
    integer :: numFlag

    write(*,*) "Input file loading..."

    !  read control parameters
    !  initial the Dx & Dy & Dz of node (2xRmax to 4xRmax).
    open (1000,FILE="systemControl.txt",STATUS='OLD',BLANK='NULL',POSITION='REWIND')

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
    read (1000,*) isFunnelWall         !  whether use Funnel Walls
    read (1000,*) isPeriodic           !  whether use Periodic function
    read (1000,*) MAX_ACC              !  maximum contact acceleration
    read (1000,*)                         
    read (1000,*) G                    !  the totle gravity
    
    numFlag = 1
    
    !  initial the contactable walls
    if (isContactWall) then
        write(*,*) 'is Contactable walls, loading...'
        read (1000,*) 
        read (1000,*) contactWallNum
        allocate (contactWallTag(contactWallNum))
        allocate (contactWallPoint(3,contactWallNum))
        allocate (contactWallVector(3,contactWallNum))
        do I = 1,contactWallNum
            contactWallTag(I) = N + numFlag
            numFlag = numFlag + 1
            read (1000,*) (contactWallPoint(K,I),K=1,3),(contactWallVector(K,I),K=1,3)
        end do
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
            funnelWallTag(I) = N + numFlag
            numFlag = numFlag + 1
            read (1000,*) (funnelWallPoint(K,I),K=1,3),(funnelWallVector(K,I),K=1,3),(funnelWallRadius(K,I),K=1,2),funnelWallLength(I)
        end do
    end if
    
    !  initial the geometry boundary of the Particle Box
    if (isPeriodic) then
        write(*,*) 'is Periodic boundary, loading...'
        read (1000,*)
        read (1000,*) (PlaSx1p(K),K=1,3),(PlaSx1v(K),K=1,3)
        read (1000,*) (PlaSx2p(K),K=1,3),(PlaSx2v(K),K=1,3)
        read (1000,*) (PlaSy1p(K),K=1,3),(PlaSy1v(K),K=1,3)
        read (1000,*) (PlaSy2p(K),K=1,3),(PlaSy2v(K),K=1,3)
        read (1000,*) LenBox
        read (1000,*) gamma
    end if
    close(1000)
    
    !  read & print the main input parameters.
    open(2000,FILE="input_points.txt",STATUS='OLD',BLANK='NULL',POSITION='REWIND')

    read (2000,*)  N                  !  the number of all particles
    read (2000,*)  m_E,m_nu           !  Youngs module; Poisson ratio
    read (2000,*)  m_mu_d,m_mu_s      !  Friciton coefficient danamic/static
    read (2000,*)  m_COR              !  Coefficient of restitution in tangential direction
    read (2000,*)  m_mu_r,m_nita_r    !  Rolling friction; Rolling damping 
    read (2000,*)  m_Beta,m_r_cut     !  Cohesion strength; Cohesion region
    read (2000,*)  m_A                !  Damping coefficient in normal direction

    !  read initial conditions from input file.
    do  I = 1,N
        read (2000,*)  Body(I),Inertia(I),(X(K,I),K=1,3),(Xdot(K,I),K=1,3),(W(K,I),K=1,3),R(I)
    end do
    X0 = X
    close(2000)
    
    !Initialize Quaternion array
    if (isQuaternion) then
        do I = 1,NMAX
            Quaternion(1,I) = 0.0D0
            Quaternion(2,I) = 0.0D0
            Quaternion(3,I) = 0.0D0
            Quaternion(4,I) = 1.0D0
        end do
    end if

    !  read properties of Saturn and Pan
    if (isPlanet) then
        write(*,*) 'is Planet system, loading...'
        open (2000,FILE="planetProperties.txt",STATUS='OLD',BLANK='NULL',POSITION='REWIND')
        read (2000,*)
        read (2000,*)  muS                 !  gravity constant of Saturn
        read (2000,*)  muP                 !  gravity constant of Pan
        read (2000,*)  omiga               !  angular velocity of Pan's initial revolution (suppose tidal locking and contrarotating, usually negative)
        read (2000,*)  (rOrig(K),K=1,3)    !  Saturn to Pan's initial position, i.e., the origin point
        close(2000)
    end if
    
    !  change position and velocity according to PBC & Shear PBC
    if (isPeriodic) then
        call periodic
    end if
    
    !  generate mesh based partition of Particles
    call meshGenerate
    
    write(*,*)  "end input"
    
    !  initial the force in case exiting overlaps at the initial-time.
    call force

    RETURN
    END