    !********************************************************************
    !     DEMbody 2.0
    !     ***********
    !
    !     Initialization of global scalars.
    !     ---------------------------------
    !
    !********************************************************************  
    SUBROUTINE zero()

    use global
    implicit none
    
    integer I,J,K

    !  Initialize parameters and set useful constants.
    Time = 0.0D0
    Tnext = 0.0D0
    Step = 0
    CheckPointTnext = 0.0D0

    write(*,*) "Hertz list initializing..."

    !  Initialize Hertz list. Linklist-Compressed.
    allocate(Head(NMAX))
    do I = 1,NMAX
        Head(I)%No = 0             !  Length of Linklist
        Head(I)%Hertz(1) = 0.0D0   !  Hertz(1)
        Head(I)%Hertz(2) = 0.0D0   !  Hertz(2)
        Head(I)%Hertz(3) = 0.0D0   !  Hertz(3)
        Head(I)%Mrot(1) = 0.0D0    !  Mrot(1)
        Head(I)%Mrot(2) = 0.0D0    !  Mrot(2)
        Head(I)%Mrot(3) = 0.0D0    !  Mrot(3)
        Head(I)%is_touching = .false.
        Head(I)%is_slipping = .false.
        Head(I)%is_rolling = .false.
        nullify(Head(I)%prev)      !  Point to Prev 
        nullify(Head(I)%next)      !  Point to Next
    end do

    return
    end