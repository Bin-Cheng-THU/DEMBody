    !********************************************************************
    !     DEMBody 4.4
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
        !  calculate linklist employing X & LatDx
        Linklist(I) = floor((X(1,I)+LatMx)/LatDx) + 1 + floor((X(2,I)+LatMy)/LatDy)*LatNx+ floor((X(3,I)+LatMz)/LatDz)*LatNx*LatNy
    end do
    !$OMP END PARALLEL DO
    
    end