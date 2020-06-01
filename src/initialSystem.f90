    !********************************************************************
    !     DEMBody 7.0
    !     ***********
    !
    !     Initialization of global scalars.
    !     ---------------------------------
    !
    !********************************************************************  
    SUBROUTINE initialSystem()

    use global
    use loadFile
    use omp_lib
    implicit none
    
    integer I,J,K,L
    integer :: wallFlag
    integer :: nRow
    integer :: idx,idy,idz
    integer :: idgx,idgy,idgz
    real(8) :: lx,ly,lz
    real(8) :: ERR
    integer :: numNode
    integer :: numGrid
    integer :: idNode
    real(8) :: ostart,oend

    !  Initialize parameters and set useful constants.
    Time = 0.0D0
    Tnext = 0.0D0
    Step = 0
    currentStep = 0
    CheckPointTnext = 0.0D0
    refreshLattice = .true.
    refreshNum = 0
    
    !  Initialize control parameters
    write(*,*) "< System file loading..."    
    open (1000,FILE="../Input/systemControl.dembody",STATUS='OLD',BLANK='NULL',POSITION='REWIND')

    read (1000,*) vsDEMBody,vsDEMBody  !  DEM version
    if (vsDEMBody .NE. VERSION) then
        write(*,*) "DEMBody Version Mismatching!"
        stop
    end if
    read (1000,*) 
    read (1000,*) LatDx,LatDy,LatDz    !  Parallel Lattice grid interval
    read (1000,*) LatNx,LatNy,LatNz    !  Parallel Lattice grid number
    read (1000,*) LatMx,LatMy,LatMz    !  Parallel Lattice grid origin
    read (1000,*) verlet               !  verlet distance
    read (1000,*)
    read (1000,*) Deltat,Tcrit         !  output time interval and termination time
    read (1000,*) CheckPointDt         !  check point time interval
    read (1000,*) Dt                   !  time step of second-step integral
    read (1000,*)                                            
    read (1000,*) isPlanet             !  whether use Planet Gravity function
    read (1000,*) isRotSystem          !  whether use Rotary System function
    read (1000,*) islocalYORP          !  whether use local YORP function
    read (1000,*) isinertialYORP       !  whether use inertial YORP function
    read (1000,*) isQuaternion         !  whether intergrate Quaternion
    read (1000,*) isContactWall        !  whether use Contactable Walls
    read (1000,*) isMovingWall         !  whether use Moving Walls
    read (1000,*) isBondedWall         !  whether use Bonded Walls
    read (1000,*) isTriMeshWall        !  whether use TriMesh Walls
    read (1000,*) isBondedTriMeshWall  !  whether use Bonded TriMesh Walls
    read (1000,*) isFunnelWall         !  whether use Funnel Walls
    read (1000,*) isPeriodic           !  whether use Periodic function
    read (1000,*) isGravBody           !  whether use Gravity Body
    read (1000,*) isSphereBody         !  whether use Sphere Body
    read (1000,*) isBiDisperse         !  whether use BiDisperse particles
    read (1000,*) isGravTriMesh        !  whether use Gravity TriMesh
    read (1000,*) MAX_ACC              !  maximum contact acceleration
    read (1000,*)                         
    read (1000,*) G(1),G(2),G(3)       !  the totle gravity    

    verlet = verlet**2                 !  verlet = (0.2Rmin)^2

    read (1000,*)
    read (1000,*)  N                  !  the number of all particles
    read (1000,*)  m_E,m_nu           !  Youngs module; Poisson ratio
    read (1000,*)  m_Ei               !  Youngs module of BiDisperse particles
    read (1000,*)  m_mu_d,m_mu_s      !  Friciton coefficient danamic/static
    read (1000,*)  m_COR              !  Coefficient of restitution
    read (1000,*)  m_Beta             !  Irregular shape
    read (1000,*)  m_c,m_r_cut        !  Cohesion strength; Cohesion region
    read (1000,*)  m_A                !  Damping coefficient in normal direction
    read (1000,*)  m_mu_r,m_nita_r    !  Rolling friction; Rolling damping
    
    wallFlag = 1
    
    !  initial the contactable walls
    if (isContactWall) then
        write(*,*) '< is Contactable walls, loading...'
        read (1000,*) 
        read (1000,*) contactWallNum
        allocate (contactWallTag(contactWallNum))
        allocate (contactWallPoint(3,contactWallNum))
        allocate (contactWallVector(3,contactWallNum))
        do I = 1,contactWallNum
            contactWallTag(I) = NMAX + wallFlag
            wallFlag = wallFlag + 1
            read (1000,*) (contactWallPoint(K,I),K=1,3),(contactWallVector(K,I),K=1,3)
        end do
    else
        read (1000,*) 
        read (1000,*)
    end if
    
    !  initial the moving walls
    if (isMovingWall) then
        write(*,*) '< is Moving walls, loading...'
        read (1000,*) 
        read (1000,*) movingWallNum
        allocate (movingWallTag(movingWallNum))
        allocate (movingWallPoint(3,movingWallNum))
        allocate (movingWallPointInit(3,movingWallNum))
        allocate (movingWallVector(3,movingWallNum))
        allocate (movingWallVelocity(3,movingWallNum))
        do I = 1,movingWallNum
            movingWallTag(I) = NMAX + wallFlag
            wallFlag = wallFlag + 1
            read (1000,*) (movingWallPointInit(K,I),K=1,3),(movingWallVector(K,I),K=1,3)
            do K = 1,3
                movingWallPoint(K,I) = movingWallPointInit(K,I)
                movingWallVelocity(K,I) = 0.0D0
            end do
        end do
        read (1000,*) movingWallTstart,movingWallTend
        read (1000,*) movingWallA,movingWallOmega
        read (1000,*) (movingWallNormal(K),K=1,3)
        !write(*,*) '< is Moving walls, loading...'
        !read (1000,*) 
        !read (1000,*) movingWallNum
        !allocate (movingWallTag(movingWallNum))
        !allocate (movingWallPoint(3,movingWallNum))
        !allocate (movingWallVector(3,movingWallNum))
        !allocate (movingWallVelocity(3,movingWallNum))
        !!  load moving walls
        !read (1000,*)
        !open (1500,File="../Input/movingWall.moving")
        !nRow = GetFileN(1500)
        !allocate (movingWallPointStore(3,movingWallNum,nRow/movingWallNum))
        !allocate (movingWallVectorStore(3,movingWallNum,nRow/movingWallNum))
        !allocate (movingWallVelocityStore(3,movingWallNum,nRow/movingWallNum))
        !do I = 1,nRow/movingWallNum
        !    do J = 1,movingWallNum
        !        read (1500,*) (movingWallPointStore(K,J,I),K=1,3),(movingWallVectorStore(K,J,I),K=1,3),(movingWallVelocityStore(K,J,I),K=1,3)
        !    end do
        !end do
        !close(1500)
        !!  initialize movingWall
        !do I = 1,movingWallNum
        !    movingWallTag(I) = NMAX + wallFlag
        !    wallFlag = wallFlag + 1
        !    do K = 1,3
        !        movingWallPoint(K,I) = movingWallPointStore(K,I,(currentStep+1))
        !        movingWallVector(K,I) = movingWallVectorStore(K,I,(currentStep+1))
        !        movingWallVelocity(K,I) = movingWallVelocityStore(K,I,(currentStep+1))
        !    end do
        !end do
    else
        read (1000,*) 
        read (1000,*)
    end if
    
    !  initial the bonded walls
    if (isBondedWall) then
        write(*,*) '< is Bonded walls, loading...'
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
            bondedWallTag(I) = NMAX + wallFlag
            wallFlag = wallFlag + 1
            read (1000,*) (bondedWallPoint(K,I),K=1,3),(bondedWallVectorN(K,I),K=1,3),(bondedWallVectorTx(K,I),K=1,3),(bondedWallVectorTy(K,I),K=1,3),bondedWallLx(I),bondedWallLy(I)
        end do
        call attitudeQ2M(bondedWallQ,bondedWallMatI,bondedWallMatB)
        !  unified input parameters
        do K = 1,3
            bondedWallWB(K) = bondedWallMatI(K,1)*bondedWallW(1) + bondedWallMatI(K,2)*bondedWallW(2) + bondedWallMatI(K,3)*bondedWallW(3)
        end do
        !  load bonded wall mesh
        read (1000,*)
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
    
    !  initial the trimesh walls
    if (isTriMeshWall) then
        write(*,*) '< is TriMesh walls, loading...'
        read (1000,*)
        read (1000,*) trimeshWallNum
        allocate (trimeshWallTag(trimeshWallNum))
        allocate (trimeshWallPoint(3,trimeshWallNum))
        allocate (trimeshWallVectorN(3,trimeshWallNum))
        allocate (trimeshWallVectorTx(3,trimeshWallNum))
        allocate (trimeshWallVectorTy(3,trimeshWallNum))
        allocate (trimeshWallLength(4,trimeshWallNum))
        !  load trimesh wall
        read (1000,*)
        open (1500,File="../Input/trimeshWall.mesh")
        do I = 1,trimeshWallNum
            trimeshWallTag(I) = NMAX + wallFlag
            wallFlag = wallFlag + 1
            read (1500,*) (trimeshWallPoint(K,I),K=1,3),(trimeshWallVectorTx(K,I),K=1,3),(trimeshWallVectorTy(K,I),K=1,3)
            !  initial normal vector
            trimeshWallVectorN(1,I) = trimeshWallVectorTx(2,I)*trimeshWallVectorTy(3,I)-trimeshWallVectorTx(3,I)*trimeshWallVectorTy(2,I)
            trimeshWallVectorN(2,I) = trimeshWallVectorTx(3,I)*trimeshWallVectorTy(1,I)-trimeshWallVectorTx(1,I)*trimeshWallVectorTy(3,I)
            trimeshWallVectorN(3,I) = trimeshWallVectorTx(1,I)*trimeshWallVectorTy(2,I)-trimeshWallVectorTx(2,I)*trimeshWallVectorTy(1,I)
            ERR = trimeshWallVectorN(1,I)*trimeshWallVectorN(1,I) + trimeshWallVectorN(2,I)*trimeshWallVectorN(2,I) + trimeshWallVectorN(3,I)*trimeshWallVectorN(3,I)
            do K = 1,3
                trimeshWallVectorN(K,I) = trimeshWallVectorN(K,I)/sqrt(ERR)
            end do
            !  initial length vector
            trimeshWallLength(1,I) = trimeshWallVectorTx(1,I)*trimeshWallVectorTx(1,I) + trimeshWallVectorTx(2,I)*trimeshWallVectorTx(2,I) + trimeshWallVectorTx(3,I)*trimeshWallVectorTx(3,I)
            trimeshWallLength(2,I) = trimeshWallVectorTx(1,I)*trimeshWallVectorTy(1,I) + trimeshWallVectorTx(2,I)*trimeshWallVectorTy(2,I) + trimeshWallVectorTx(3,I)*trimeshWallVectorTy(3,I)
            trimeshWallLength(3,I) = trimeshWallVectorTy(1,I)*trimeshWallVectorTy(1,I) + trimeshWallVectorTy(2,I)*trimeshWallVectorTy(2,I) + trimeshWallVectorTy(3,I)*trimeshWallVectorTy(3,I)
            trimeshWallLength(4,I) = trimeshWallLength(2,I)*trimeshWallLength(2,I) - trimeshWallLength(1,I)*trimeshWallLength(3,I)
        end do
        close(1500)
        !  Distribute into Lattice
        call initialTriMesh
    else
        read (1000,*)
        read (1000,*)
    end if

    !  initial the bonded trimesh walls
    if (isBondedTriMeshWall) then
        write(*,*) '< is Bonded TriMesh walls, loading...'
        read (1000,*)
        read (1000,*) (bondedTriMeshWallX(K),K=1,3),(bondedTriMeshWallXdot(K),K=1,3),(bondedTriMeshWallW(K),K=1,3)
        read (1000,*) (bondedTriMeshWallQ(K),K=1,4)
        read (1000,*) bondedTriMeshWallBody,(bondedTriMeshWallInertia(K),K=1,3)
        read (1000,*) bondedTriMeshWallNum
        allocate (bondedTriMeshWallTag(bondedTriMeshWallNum))
        allocate (bondedTriMeshWallPoint(3,bondedTriMeshWallNum))
        allocate (bondedTriMeshWallVectorN(3,bondedTriMeshWallNum))
        allocate (bondedTriMeshWallVectorTx(3,bondedTriMeshWallNum))
        allocate (bondedTriMeshWallVectorTy(3,bondedTriMeshWallNum))
        allocate (bondedTriMeshWallLength(4,bondedTriMeshWallNum))
        !  load trimesh wall
        read (1000,*)
        open (1500,File="../Input/bondedTriMeshWall.mesh")
        do I = 1,bondedTriMeshWallNum
            bondedTriMeshWallTag(I) = NMAX + wallFlag
            wallFlag = wallFlag + 1
            read (1500,*) (bondedTriMeshWallPoint(K,I),K=1,3),(bondedTriMeshWallVectorTx(K,I),K=1,3),(bondedTriMeshWallVectorTy(K,I),K=1,3)
            !  initial normal vector
            bondedTriMeshWallVectorN(1,I) = bondedTriMeshWallVectorTx(2,I)*bondedTriMeshWallVectorTy(3,I)-bondedTriMeshWallVectorTx(3,I)*bondedTriMeshWallVectorTy(2,I)
            bondedTriMeshWallVectorN(2,I) = bondedTriMeshWallVectorTx(3,I)*bondedTriMeshWallVectorTy(1,I)-bondedTriMeshWallVectorTx(1,I)*bondedTriMeshWallVectorTy(3,I)
            bondedTriMeshWallVectorN(3,I) = bondedTriMeshWallVectorTx(1,I)*bondedTriMeshWallVectorTy(2,I)-bondedTriMeshWallVectorTx(2,I)*bondedTriMeshWallVectorTy(1,I)
            ERR = bondedTriMeshWallVectorN(1,I)*bondedTriMeshWallVectorN(1,I) + bondedTriMeshWallVectorN(2,I)*bondedTriMeshWallVectorN(2,I) + bondedTriMeshWallVectorN(3,I)*bondedTriMeshWallVectorN(3,I)
            do K = 1,3
                bondedTriMeshWallVectorN(K,I) = bondedTriMeshWallVectorN(K,I)/sqrt(ERR)
            end do
            !  initial length vector
            bondedTriMeshWallLength(1,I) = bondedTriMeshWallVectorTx(1,I)*bondedTriMeshWallVectorTx(1,I) + bondedTriMeshWallVectorTx(2,I)*bondedTriMeshWallVectorTx(2,I) + bondedTriMeshWallVectorTx(3,I)*bondedTriMeshWallVectorTx(3,I)
            bondedTriMeshWallLength(2,I) = bondedTriMeshWallVectorTx(1,I)*bondedTriMeshWallVectorTy(1,I) + bondedTriMeshWallVectorTx(2,I)*bondedTriMeshWallVectorTy(2,I) + bondedTriMeshWallVectorTx(3,I)*bondedTriMeshWallVectorTy(3,I)
            bondedTriMeshWallLength(3,I) = bondedTriMeshWallVectorTy(1,I)*bondedTriMeshWallVectorTy(1,I) + bondedTriMeshWallVectorTy(2,I)*bondedTriMeshWallVectorTy(2,I) + bondedTriMeshWallVectorTy(3,I)*bondedTriMeshWallVectorTy(3,I)
            bondedTriMeshWallLength(4,I) = bondedTriMeshWallLength(2,I)*bondedTriMeshWallLength(2,I) - bondedTriMeshWallLength(1,I)*bondedTriMeshWallLength(3,I)
        end do
        close(1500)
        !  Distribute into Lattice
        call initialBondedTriMesh
    else
        read (1000,*)
        read (1000,*)
    end if
    
    !  initial the funnel walls
    if (isFunnelWall) then
        write(*,*) '< is Funnel walls, loading...'
        read (1000,*) 
        read (1000,*) funnelWallNum
        allocate (funnelWallTag(funnelWallNum))
        allocate (funnelWallPoint(3,funnelWallNum))
        allocate (funnelWallVector(3,funnelWallNum))
        allocate (funnelWallRadius(2,funnelWallNum))
        allocate (funnelWallLength(funnelWallNum))
        do I = 1,funnelWallNum
            funnelWallTag(I) = NMAX + wallFlag
            wallFlag = wallFlag + 1
            read (1000,*) (funnelWallPoint(K,I),K=1,3),(funnelWallVector(K,I),K=1,3),(funnelWallRadius(K,I),K=1,2),funnelWallLength(I)
        end do
    else
        read (1000,*)
        read (1000,*)
    end if
    
    !  initial the geometry boundary of the Particle Box
    if (isPeriodic) then
        write(*,*) '< is Periodic boundary, loading...'
        read (1000,*)
        read (1000,*) (PlaSx1p(K),K=1,3),(PlaSx1v(K),K=1,3)
        read (1000,*) (PlaSx2p(K),K=1,3),(PlaSx2v(K),K=1,3)
        read (1000,*) (PlaSy1p(K),K=1,3),(PlaSy1v(K),K=1,3)
        read (1000,*) (PlaSy2p(K),K=1,3),(PlaSy2v(K),K=1,3)
        read (1000,*) LenBoxX
        read (1000,*) LenBoxY
        read (1000,*) gamma
        LenBoxX = abs(PlaSx1p(1) - PlaSx2p(1))
        LenBoxY = abs(PlaSy1p(2) - PlaSy2p(2))
    else
        read (1000,*)
        read (1000,*)
    end if
    
    !  initial the gravity body
    if (isGravBody) then
        write(*,*) '< is Gravity Body, loading...'
        read (1000,*)
        gravBodyTag = NMAX + wallFlag
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
    
    !  initial the sphere body
    if (isSphereBody) then
        write(*,*) '< is Sphere Body, loading...'
        read (1000,*)
        read (1000,*) sphereBodyNum
        allocate (sphereBodyTag(sphereBodyNum))
        allocate (sphereBodyX(3,sphereBodyNum))
        allocate (sphereBodyXdot(3,sphereBodyNum))
        allocate (sphereBodyW(3,sphereBodyNum))
        allocate (sphereBodyQ(4,sphereBodyNum))
        allocate (sphereBodyBody(sphereBodyNum))
        allocate (sphereBodyR(sphereBodyNum))
        allocate (sphereBodyInertia(sphereBodyNum))
        allocate (sphereBodyF(3,sphereBodyNum))
        allocate (sphereBodyFM(3,sphereBodyNum))
        do I = 1,sphereBodyNum
            sphereBodyTag(I) = NMAX + wallFlag
            wallFlag = wallFlag + 1
            read (1000,*) sphereBodyBody(I),sphereBodyInertia(I),(sphereBodyX(K,I),K=1,3),(sphereBodyXdot(K,I),K=1,3),(sphereBodyW(K,I),K=1,3),sphereBodyR(I)
            sphereBodyQ(1,I) = 0.0D0
            sphereBodyQ(2,I) = 0.0D0
            sphereBodyQ(3,I) = 0.0D0
            sphereBodyQ(4,I) = 1.0D0
        end do
    else
        read (1000,*)
        read (1000,*)
    end if

    !  initial the biDisperse body
    if (isBiDisperse) then
        write(*,*) '< is Bi-Disperse, loading...'
        read (1000,*)
        read (1000,*) biDisperseNum
        !read (1000,*) biDisperseScale  !  biDisperseR/LatDx == biDisperseR/Rmax/2.5

        read (1000,*) verletBiDisperse  !  0.75Rmax
        verletBiDisperse = verletBiDisperse**2
       
        allocate (biDisperseTag(biDisperseNum))
        allocate (biDisperseX(3,biDisperseNum))
        allocate (biDisperseXdot(3,biDisperseNum))
        allocate (biDisperseW(3,biDisperseNum))
        allocate (biDisperseQ(4,biDisperseNum))
        allocate (biDisperseBody(biDisperseNum))
        allocate (biDisperseR(biDisperseNum))
        allocate (biDisperseInertia(biDisperseNum))
        allocate (biDisperseF(3,biDisperseNum))
        allocate (biDisperseFM(3,biDisperseNum))
        allocate (biDisperseXT(3,biDisperseNum))
        allocate (biDisperseEnergy(biDisperseNum))
        allocate (biDisperseHeat(biDisperseNum))
        !  load bi-disperse particles data
        read (1000,*)
        open (1500,File="../Input/largeParticles.bidisperse")
        do I = 1,biDisperseNum
            biDisperseTag(I) = NMAX + wallFlag
            wallFlag = wallFlag + 1
            read (1500,*) biDisperseBody(I),biDisperseInertia(I),(biDisperseX(K,I),K=1,3),(biDisperseXdot(K,I),K=1,3),(biDisperseW(K,I),K=1,3),biDisperseR(I)
            biDisperseQ(1,I) = 0.0D0
            biDisperseQ(2,I) = 0.0D0
            biDisperseQ(3,I) = 0.0D0
            biDisperseQ(4,I) = 1.0D0
        end do
        close(1500)
        call initialBiDisperse
    else
        read (1000,*)
        read (1000,*)
    end if
    
    !  initial the properties of Saturn and Pan
    if (isPlanet) then
        write(*,*) '< is Planet system, loading'
        read (1000,*)
        read (1000,*) muS, muP            !  Gravity constant of Saturn and Pan
        read (1000,*) omega               !  Angular velocity of Pan's initial revolution (suppose tidal locking and anticlockwise, usually negative)
        read (1000,*) (rOrig(K),K=1,3)    !  Saturn to Pan's initial position, i.e., the origin point
        omega = int(omega/abs(omega))*sqrt(abs(muS/rOrig(2)/rOrig(2)/rOrig(2)))
    else
        read (1000,*)
        read (1000,*)
    end if
        
    !  initial the rotary system
    if (isRotSystem) then
        write(*,*) '< is Rotary System, loading...'
        read (1000,*)
        read (1000,*) sysOmega
    else
        read (1000,*)
        read (1000,*)
    end if

    !  initial the localYORP system
    if (islocalYORP) then
        write(*,*) '< is YORP System in local frame, loading...'
        read (1000,*)
        read (1000,*) localYORPTstart,localYORPTend
        read (1000,*) localYORPOmega
        read (1000,*) localYORPIncrement,localYORPDt
        localYORPTnext = localYORPTstart
        read (1000,*) localYORPdensity
        read (1000,*) localYORPradius
        localYORPmass = 4.0D0/3.0D0*PI*localYORPradius**3*localYORPdensity
        localYORPcenter = 0.0D0
    else
        read (1000,*)
        read (1000,*)
    end if

    !  initial the inertialYORP system
    if (isinertialYORP) then
        write(*,*) '< is YORP System in inertial frame, loading...'
        read (1000,*)
        read (1000,*) inertialYORPTstart,inertialYORPTend
        read (1000,*) inertialYORPOmega
        read (1000,*) inertialYORPIncrement,inertialYORPDt
        inertialYORPTnext = inertialYORPTstart
        read (1000,*) inertialYORPdensity
        read (1000,*) inertialYORPradius
        read (1000,*) inertialYORPAllmass
        inertialYORPmass = 4.0D0/3.0D0*PI*inertialYORPradius**3*inertialYORPdensity
        inertialYORPcenter = 0.0D0
    else
        read (1000,*)
        read (1000,*)
    end if

    !  initial the grav trimesh
    if (isGravTriMesh) then
        write(*,*) '< is Grav TriMesh, loading...'
        read (1000,*)
        read (1000,*) TriLatDx,TriLatDy,TriLatDz
        read (1000,*) TriLatMx,TriLatMy,TriLatMz
        read (1000,*) TriLatNx,TriLatNy,TriLatNz
        read (1000,*)
        open (2000,FILE="../Input/gravTriMesh.force",STATUS='OLD',BLANK='NULL',POSITION='REWIND')
        numNode = (TriLatNx+1)*(TriLatNy+1)*(TriLatNz+1)
        numGrid = TriLatNx*TriLatNy*TriLatNz
        allocate(GravTriMeshNode(3,numNode))
        do I = 1,numNode
            read (2000,*) (GravTriMeshNode(K,I),K=1,3)
        end do
        close(2000)
        
        allocate(GravTriMeshGrid(numGrid))
        do I = 1,numGrid
            idz = int((I-1)/(TriLatNx*TriLatNy))+1
            idy = int((I-(idz-1)*(TriLatNx*TriLatNy)-1)/TriLatNx)+1
            idx = I-(idz-1)*(TriLatNx*TriLatNy)-(idy-1)*TriLatNx
            lx = dble(idx-1)*TriLatDx - TriLatMx
            ly = dble(idy-1)*TriLatDy - TriLatMy
            lz = dble(idz-1)*TriLatDz - TriLatMz            
            GravTriMeshGrid(I)%PositionCenter(1) = lx + 0.5D0*TriLatDx
            GravTriMeshGrid(I)%PositionCenter(2) = ly + 0.5D0*TriLatDy
            GravTriMeshGrid(I)%PositionCenter(3) = lz + 0.5D0*TriLatDz
            idNode = idx + (idy-1)*(TriLatNx+1) + (idz-1)*(TriLatNx+1)*(TriLatNy+1)
            GravTriMeshGrid(I)%NodeID(1) = idNode
            GravTriMeshGrid(I)%NodeID(2) = idNode + 1
            GravTriMeshGrid(I)%NodeID(3) = idNode + 2 + TriLatNx
            GravTriMeshGrid(I)%NodeID(4) = idNode + 1 + TriLatNx
            GravTriMeshGrid(I)%NodeID(5) = idNode + (TriLatNx+1)*(TriLatNy+1)
            GravTriMeshGrid(I)%NodeID(6) = idNode + 1 + (TriLatNx+1)*(TriLatNy+1)
            GravTriMeshGrid(I)%NodeID(7) = idNode + 2 + TriLatNx + (TriLatNx+1)*(TriLatNy+1)
            GravTriMeshGrid(I)%NodeID(8) = idNode + 1 + TriLatNx + (TriLatNx+1)*(TriLatNy+1)
        end do
    else
        read (1000,*)
        read (1000,*)
    end if
    close(1000)

    !  initial DEM Lattice and Periodic DEM Lattice
    call initialLattice
    
    return
    end