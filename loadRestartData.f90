    !********************************************************************
    !     DEMbody 2.0
    !     ***********
    !
    !     Load restart data including Hertz contact history.
    !     --------------------------------------------------
    !
    !********************************************************************
    subroutine loadRestartData()

    use loadFile
    Implicit None

    Character( Len = 2048 ) :: cLine
    integer :: nRow , nCol , I
    integer K  
    logical alive

    write(*,*) "RestartData loading..."

    INQUIRE( File = 'RestartData.txt', EXIST = alive)
    IF (alive) then
        Open( 12 , File = 'RestartData.txt' )
        nRow = GetFileN( 12 )
        write( * , * ) 'RestartData ',nRow,' rows!'

        Do I = 1 , nRow
            read( 12 , '(a2048)' ) cLine
            BACKSPACE(12)
            nCol = GetDataN( cLine )
            call Load(12, nCol)
        End Do

        close( 12 )

    else
        write(*,*) "No RestartData exist."
        write(*,*) "Go to calculating directly."
    end if

    End subroutine