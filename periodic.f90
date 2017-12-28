    !********************************************************************
    !     DEMbody 2.0
    !     ***********
    !
    !     Change position and velocity when crossing boundary.
    !     --------------------------
    !
    !     @Periodic boundary
    !     @Shear periodic boundary
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
            X(1,I) = X(1,I) - LenBox
        end if
        !Wallx2
        if (X(1,I).LT.(PlaSx2p(1))) then
            X(1,I) = X(1,I) + LenBox
        end if
        !Wally1
        if (X(2,I).GT.(PlaSy1p(2))) then
            X(2,I) = X(2,I) - LenBox
            X(1,I) = X(1,I) + 0.5*LenBox*gamma*time
            X(1,I) = X(1,I) - ((X(1,I)+0.5*LenBox) - mod((X(1,I)+0.5*LenBox),LenBox))
            Xdot(1,I) = Xdot(1,I) + LenBox*gamma
        end if
        !Wally2
        if (X(2,I).LT.(PlaSy2p(2))) then
            X(2,I) = X(2,I) + LenBox
            X(1,I) = X(1,I) - 0.5*LenBox*gamma*time
            X(1,I) = X(1,I) - ((X(1,I)-0.5*LenBox) - mod((X(1,I)-0.5*LenBox),LenBox))
            Xdot(1,I) = Xdot(1,I) - LenBox*gamma
        end if

        if ((PlaSx1p(1)-X(1,I)-2D0*R(I)).LE.0.0) then
            Tag2(I) = 1
        else
            Tag2(I) = 0
        end if

        if ((X(1,I)-PlaSx2p(1)-2D0*R(I)).LE.0.0) then
            Tag1(I) = 1
        else
            Tag1(I) = 0
        end if

        if ((PlaSy1p(2)-X(2,I)-2D0*R(I)).LE.0.0) then
            Tag4(I) = 1
        else
            Tag4(I) = 0
        end if

        if ((X(2,I)-PlaSy2p(2)-2D0*R(I)).LE.0.0) then
            Tag3(I) = 1
        else
            Tag3(I) = 0
        end if

    end do
    !$OMP END PARALLEL DO

    end