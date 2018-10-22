    !********************************************************************
    !     DEMBody 4.4
    !     ***********
    !
    !     Force for all particles.
    !     --------------------------
    !      
    !     @Using forceParticleLattice or forceLattice when enable Lattice
    !     @Using forceTraverse when enable Traverse
    !     @Using forceMirror when enable Periodic
    !     @Using forceContactWalls when enable ContactWall
    !     @Using forceBondedWalls when enable BondedWall
    !     @Using forceTriMeshWalls when enable TriMeshWall
    !     @Using forceFunnelWalls when enable FunnelWall
    !     @Using forceGravBody when enable GravBody
    !     @force to acceleration
    !     @Using Planet force when enable Planet
    !     @Using forceRotSystem when enbale RotSystem
    !     @Using forceGravTriMesh when enable GravTriMesh
    !
    !********************************************************************
    subroutine force()

    use global
    use omp_lib
    implicit none
    
    real(8) :: ostart,oend
    integer :: I,J,K
    real(8) :: accMag
    character(30) :: FileNameHead
    character(30) :: FileNameForce

    !################         Part 1: Particle forces          ###################
#ifdef LatticeSearch
#ifdef ParticleLattice
    call forceParticleLattice
#else
    call forceLattice
#endif
#elif TraverseSearch
    call forceTraverse
#endif

    !################         Part 2: Mirrored Particle forces          ###################
    !ostart = omp_get_wtime()
    !  calculate mirrored force if using PBC & Shear PBC
    if (isPeriodic) then
        call forceMirror
    end if
    !oend = omp_get_wtime()
    !write(*,*) "Mirror",(oend-ostart)
    
    !################         Part 3: Contact wall forces          ###################
    !ostart = omp_get_wtime()
    !  calculate force of contactable walls if using contactable walls
    if (isContactWall) then
        call forceContactWalls
    end if
    !oend = omp_get_wtime()
    !write(*,*) "Contact walls", (oend-ostart)
    
    !################         Part 4: Bonded wall forces          ###################
    !  calculate force of bonded walls if using bonded walls
    !ostart = omp_get_wtime()
    !  calculate force of bonded walls if using bonded walls
    if (isBondedWall) then                     
        call forceBondedWalls 
    end if
    !oend = omp_get_wtime()
    !write(*,*) "Bonded walls",(oend-ostart)
    
    !################         Part 5: Trimesh wall forces          ###################
    !ostart = omp_get_wtime()
    !  calculate force of trimesh walls if using trimesh walls
    if (isTriMeshWall) then                     
        call forceTriMeshWalls
    end if
    !oend = omp_get_wtime()
    !write(*,*) "Tri meshes",(oend-ostart)
    
    !################         Part 6: Funnel wall forces          ###################
    !ostart = omp_get_wtime()
    !  calculate force of funnel walls if using funnel walls
    if (isFunnelWall) then
        call forceFunnelWalls
    end if
    !oend = omp_get_wtime()
    !write(*,*) "Funnel walls",(oend-ostart)
    
    !################         Part 7: GravBody forces          ###################
    !ostart = omp_get_wtime()
    !  calculate force of gravity body if using GravBody
    if (isGravBody) then
        call forceGravBody
    end if
    !oend = omp_get_wtime()
    !write(*,*) "Grav bodies",(oend-ostart)
    
    !################         Part 8: Force to Acceleration          ###################
    !ostart = omp_get_wtime()
    !$OMP PARALLEL DO PRIVATE(I,K,accMag)
    do I = 1,N
        do K = 1,3
            F(K,I) = F(K,I)/Body(I)
            FM(K,I) = FM(K,I)/Inertia(I)
        end do
        accMag = sqrt(F(1,I)**2 + F(2,I)**2+F(3,I)**2)
        if (accMag .GE. MAX_ACC) then
            write(*,*) 'exceed Max Acc',I,accMag
            do K = 1,3
                F(K,I) = F(K,I)/accMag*MAX_ACC
            end do
        end if
        F(3,I) = F(3,I) - G
    end do
    !$OMP END PARALLEL DO
    !oend = omp_get_wtime()
    !write(*,*) 'accelerate', (oend-ostart)    

    !################         Part 9: Planet forces          ###################
    !ostart = omp_get_wtime()
    !  calculate Planet force if in Planet system
    if (isPlanet) then
        call Planet
    end if
    !oend = omp_get_wtime()
    !write(*,*) 'Planet', (oend-ostart)  
    
    !################         Part 10: Rotsystem forces          ###################
    !ostart = omp_get_wtime()
    !  calculate Noninertial force if in Rotary system
    if (isRotSystem) then
        call forceRotSystem
    end if
    !oend = omp_get_wtime()
    !write(*,*) 'Rot system', (oend-ostart)  
    
    !################         Part 11: GravTriMesh forces          ###################
    !ostart = omp_get_wtime()
    !  calculate force of gravity body if using GravBody
    if (isGravTriMesh) then
        call forceGravTriMesh
    end if
    !oend = omp_get_wtime()
    !write(*,*) "Grav bodies",(oend-ostart)
    
    !write(FileNameForce,'(F10.5)') Time
    !FileNameForce = trim(FileNameForce)//'Force.txt'
    !open(10000,FILE=FileNameForce)
    !do I = 1,N
    !    write(10000,'(3F30.15)') (F(K,I),K=1,3)
    !end do
    !close(10000)
    !close(123)
    !pause
    !
    !stop
    
    return
    end
