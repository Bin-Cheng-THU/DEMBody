    !********************************************************************
    !     DEMBody 6.0
    !     ***********
    !
    !     Change position and velocity when crossing boundary.
    !     --------------------------
    !
    !     @Periodic boundary
    !     @Shear periodic boundary
    !     @Can not solve shear periodic boundary in current version (v6.0)
    !
    !********************************************************************
    subroutine periodic()

    use global
    implicit none

    integer :: I
    
    !$OMP PARALLEL DO PRIVATE(I)
    do I = 1,N
        !Wallx1
        if (X(1,I).GT.(PlaSx1p(1))) then
            X(1,I) = X(1,I) - LenBoxX
        end if
        !Wallx2
        if (X(1,I).LT.(PlaSx2p(1))) then
            X(1,I) = X(1,I) + LenBoxX
        end if
        !Wally1
        if (X(2,I).GT.(PlaSy1p(2))) then
            X(2,I) = X(2,I) - LenBoxY
        end if
        !Wally2
        if (X(2,I).LT.(PlaSy2p(2))) then
            X(2,I) = X(2,I) + LenBoxY
        end if
    end do
    !$OMP END PARALLEL DO

    end