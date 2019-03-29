    !********************************************************************
    !     DEMBody 5.1
    !     ***********
    !
    !     Initialization of biDisperse.
    !     ---------------------------------
    !
    !********************************************************************
    SUBROUTINE initialBiDisperse()

    use global
    implicit none

    integer I,J,K,L

    !  initialize biDisperse
    LatNum = LatNx*LatNy*LatNz 
    allocate(biDisperseDEM(LatNum))
    allocate(biDisperseDEMtail(LatNum))
    do I = 1,LatNum
        biDisperseDEM(I)%No = 0
        nullify(biDisperseDEM(I)%next)
    end do
    
    !  initialize  HeadBiDisperse
    allocate(HeadBiDisperse(biDisperseNum))
    do I = 1,biDisperseNum
        HeadBiDisperse(I)%No = 0                    !  Length of Linklist
        !HeadBiDisperse(I)%recordTime = 0.0D0        !  record time
        HeadBiDisperse(I)%Hertz(1) = 0.0D0          !  Hertz(1)
        HeadBiDisperse(I)%Hertz(2) = 0.0D0          !  Hertz(2)
        HeadBiDisperse(I)%Hertz(3) = 0.0D0          !  Hertz(3)
        HeadBiDisperse(I)%Mrot(1) = 0.0D0           !  Mrot(1)
        HeadBiDisperse(I)%Mrot(2) = 0.0D0           !  Mrot(2)
        HeadBiDisperse(I)%Mrot(3) = 0.0D0           !  Mrot(3)
        HeadBiDisperse(I)%Mtwist(1) = 0.0D0         !  Mtwist(1)
        HeadBiDisperse(I)%Mtwist(2) = 0.0D0         !  Mtwist(2)
        HeadBiDisperse(I)%Mtwist(3) = 0.0D0         !  Mtwist(3)
        HeadBiDisperse(I)%is_touching = .false.     !  Whether touch or not
        HeadBiDisperse(I)%is_slipping = .false.     !  Whether slip or not
        HeadBiDisperse(I)%is_rolling = .false.      !  Whether roll or not
        HeadBiDisperse(I)%is_twisting = .false.     !  Whether twist or not
        nullify(HeadBiDisperse(I)%prev)             !  Point to Prev 
        nullify(HeadBiDisperse(I)%next)             !  Point to Next
    end do
    
    end