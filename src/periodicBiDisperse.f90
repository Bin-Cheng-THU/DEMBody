    !********************************************************************
    !     DEMBody 5.2
    !     ***********
    !
    !     Change position and velocity when crossing boundary.
    !     --------------------------
    !
    !     @Periodic boundary
    !     @Shear periodic boundary
    !
    !********************************************************************
    subroutine periodicBiDisperse()

    use global
    implicit none

    integer :: I
    
    !$OMP PARALLEL DO PRIVATE(I)
    do I = 1,biDisperseNum
        !Wallx1
        if (biDisperseX(1,I).GT.(PlaSx1p(1))) then
            biDisperseX(1,I) = biDisperseX(1,I) - LenBoxX
            refreshLattice = .true.
        end if
        !Wallx2
        if (biDisperseX(1,I).LT.(PlaSx2p(1))) then
            biDisperseX(1,I) = biDisperseX(1,I) + LenBoxX
            refreshLattice = .true.
        end if
        !Wally1
        if (biDisperseX(2,I).GT.(PlaSy1p(2))) then
            biDisperseX(2,I) = biDisperseX(2,I) - LenBoxY
            biDisperseX(1,I) = biDisperseX(1,I) + LenBoxY*gamma*time
            biDisperseX(1,I) = biDisperseX(1,I) - ((biDisperseX(1,I)+0.5D0*LenBoxX) - MODULO((biDisperseX(1,I)+0.5D0*LenBoxX),LenBoxX))
            biDisperseXdot(1,I) = biDisperseXdot(1,I) + LenBoxY*gamma
            refreshLattice = .true.
        end if
        !Wally2
        if (biDisperseX(2,I).LT.(PlaSy2p(2))) then
            biDisperseX(2,I) = biDisperseX(2,I) + LenBoxY
            biDisperseX(1,I) = biDisperseX(1,I) - LenBoxY*gamma*time
            biDisperseX(1,I) = biDisperseX(1,I) - ((biDisperseX(1,I)+0.5D0*LenBoxX) - MODULO((biDisperseX(1,I)+0.5D0*LenBoxX),LenBoxX))
            biDisperseXdot(1,I) = biDisperseXdot(1,I) - LenBoxY*gamma
            refreshLattice = .true.
        end if

        !  Periodic Tag for biDisperse particles
        if ((PlaSx1p(1)-biDisperseX(1,I)-LatDx*(biDisperseScale+1)).LE.0.0D0) then
            BiLTag2(I) = 1
        else
            BiLTag2(I) = 0
        end if
        if ((biDisperseX(1,I)-PlaSx2p(1)-LatDx*(biDisperseScale+1)).LE.0.0D0) then
            BiLTag1(I) = 1
        else
            BiLTag1(I) = 0
        end if
        if ((PlaSy1p(2)-biDisperseX(2,I)-LatDy*(biDisperseScale+1)).LE.0.0D0) then
            BiLTag4(I) = 1
        else
            BiLTag4(I) = 0
        end if
        if ((biDisperseX(2,I)-PlaSy2p(2)-LatDy*(biDisperseScale+1)).LE.0.0D0) then
            BiLTag3(I) = 1
        else
            BiLTag3(I) = 0
        end if
        
        !  Periodic Tag for inner biDisperse particles
        if ((PlaSx1p(1)-biDisperseX(1,I)-LatDx*(2*biDisperseScale+1)).LE.0.0D0) then
            BiTag2(I) = 1
        else
            BiTag2(I) = 0
        end if
        if ((biDisperseX(1,I)-PlaSx2p(1)-LatDx*(2*biDisperseScale+1)).LE.0.0D0) then
            BiTag1(I) = 1
        else
            BiTag1(I) = 0
        end if
        if ((PlaSy1p(2)-biDisperseX(2,I)-LatDy*(2*biDisperseScale+1)).LE.0.0D0) then
            BiTag4(I) = 1
        else
            BiTag4(I) = 0
        end if
        if ((biDisperseX(2,I)-PlaSy2p(2)-LatDy*(2*biDisperseScale+1)).LE.0.0D0) then
            BiTag3(I) = 1
        else
            BiTag3(I) = 0
        end if
    end do
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO PRIVATE(I)
    do I = 1,N
        !  Periodic Tag for particles when checking contact with bidisperse particles
        if ((PlaSx1p(1)-X(1,I)-LatDx*(biDisperseScale+1)).LE.0.0D0) then
            BiSTag2(I) = 1
        else
            BiSTag2(I) = 0
        end if

        if ((X(1,I)-PlaSx2p(1)-LatDx*(biDisperseScale+1)).LE.0.0D0) then
            BiSTag1(I) = 1
        else
            BiSTag1(I) = 0
        end if

        if ((PlaSy1p(2)-X(2,I)-LatDy*(biDisperseScale+1)).LE.0.0D0) then
            BiSTag4(I) = 1
        else
            BiSTag4(I) = 0
        end if

        if ((X(2,I)-PlaSy2p(2)-LatDy*(biDisperseScale+1)).LE.0.0D0) then
            BiSTag3(I) = 1
        else
            BiSTag3(I) = 0
        end if
    end do
    !$OMP END PARALLEL DO

    end