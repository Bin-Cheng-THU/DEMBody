    !********************************************************************
    !     DEMBody 4.6
    !     ***********
    !
    !     Parameter input.
    !     ----------------
    !
    !********************************************************************
    subroutine loadData()

    use global
    implicit none
    
    integer :: I,K

    write(*,*) "< Input file loading..."
    
    !  read & print the main input parameters.
    open (2000,FILE="../Input/input_points.txt",STATUS='OLD',BLANK='NULL',POSITION='REWIND')

    read (2000,*)  N                  !  the number of all particles
    read (2000,*)  m_E,m_nu           !  Youngs module; Poisson ratio
    read (2000,*)  m_mu_d,m_mu_s      !  Friciton coefficient danamic/static
    read (2000,*)  m_COR              !  Coefficient of restitution
    read (2000,*)  m_Beta             !  Irregular shape
    read (2000,*)  m_c,m_r_cut        !  Cohesion strength; Cohesion region
    read (2000,*)  m_A                !  Damping coefficient in normal direction
    read (2000,*)  m_mu_r,m_nita_r    !  Rolling friction; Rolling damping

    !  read initial conditions from input file.
    do  I = 1,N
        read (2000,*)  Body(I),Inertia(I),(X(K,I),K=1,3),(Xdot(K,I),K=1,3),(W(K,I),K=1,3),R(I)
    end do
    X0 = X
    close(2000)
    
    !  initialize Quaternion array
    if (isQuaternion) then
        do I = 1,NMAX
            Quaternion(1,I) = 0.0D0
            Quaternion(2,I) = 0.0D0
            Quaternion(3,I) = 0.0D0
            Quaternion(4,I) = 1.0D0
        end do
    end if
    
    !  change position and velocity according to PBC & Shear PBC
    if (isPeriodic) then
        call periodic
    end if

#ifdef LatticeSearch
    !  generate mesh based on partition of Particles
    call meshGenerate
    
    !  generate DEM Parallel Latice based on partition of Particles
    call latticeGenerate
#endif

    !  initial the force in case exiting overlaps at the initial-time.
    call force

    RETURN
    END