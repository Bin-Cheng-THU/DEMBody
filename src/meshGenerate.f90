    !********************************************************************
    !     DEMBody 4.3
    !     ***********
    !
    !     Generate mesh based on partition of Particles.
    !     -------------------------------------------
    !
    !********************************************************************
    subroutine meshGenerate()

    use global
    implicit none

    integer :: I

    !$OMP PARALLEL DO PRIVATE(I)
    do I = 1,N
        !  calculate linklist employing X & Dx
        Linklist(I) = floor((X(1,I)+Mx)/Dx) + Nx + floor((X(2,I)+My)/Dy)*Ny+ floor((X(3,I)+Mz)/Dz)*Nz
    end do
    !$OMP END PARALLEL DO
    
    end