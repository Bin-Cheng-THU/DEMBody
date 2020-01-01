    !********************************************************************
    !     DEMBody 6.0
    !     ***********
    !
    !     Parameter input.
    !     ----------------
    !
    !********************************************************************
    subroutine loadData()

    use global
    implicit none
    
    integer :: I,J,K

    write(*,*) "< Input file loading..."
    
    !  read & print the main input parameters.
    open (2000,FILE="../Input/input_points.txt",STATUS='OLD',BLANK='NULL',POSITION='REWIND')

    !  read initial conditions from input file.
    do I = 1,N
        read (2000,*)  Body(I),Inertia(I),(bondParticles(K,I),K=1,3),R(I)
    end do
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
    
    !  read & print the main input parameters.
    open (3000,FILE="../Input/input_assembly.txt",STATUS='OLD',BLANK='NULL',POSITION='REWIND')
    read (3000,*)  bondN              !  the number of all assemblies
    call initialAssembly
    close(3000)
    
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