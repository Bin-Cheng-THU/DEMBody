    !********************************************************************
    !     DEMbody 2.0
    !     ***********
    !
    !     Generate mesh based partition of Particles.
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
        Linklist(I) = int((X(1,I)+Mx)/Dx) + Nx + int((X(2,I)+My)/Dy)*Ny+ int((X(3,I)+Mz)/Dz)*Nz
    end do
    !$OMP END PARALLEL DO
    
    end