    !********************************************************************
    !     DEMBody 4.6
    !     ***********
    !
    !     Load file modules.
    !     --------------------------
    !   
    !     @GetDataN()      Get columns of rows
    !     @Load()          Load data to Hertz linklist
    !     @GetFileN()      Get rows of files
    !
    !********************************************************************
    module loadFile
    use global
    implicit None

    contains 

    Integer Function GetDataN( cStr )
    Character( Len = * ) , Intent( IN ) :: cStr
    Integer :: i
    Logical :: bIsSeparator , bIsQuote
    GetDataN = 0
    bIsSeparator = .TRUE.
    bIsQuote = .FALSE.
    Do i = 1 , Len_Trim( cStr )
        Select Case( cStr(i:i) )
        Case( '"' , "'" ) !// quotation mark
            If ( .Not.bIsQuote ) GetDataN = GetDataN + 1  !//if not in quotation
            bIsQuote = .Not.bIsQuote !// begaining or ending of quotation
            bIsSeparator = .FALSE.
        Case( " " , "," , char(9) ) !// comma of blank
            If ( .Not.bIsQuote ) then  !// if not in quotation
                bIsSeparator = .TRUE.
            End If
            Case Default      
            If ( bIsSeparator ) then
                GetDataN = GetDataN + 1
            End If
            bIsSeparator = .FALSE.
        End Select
    End Do
    End Function GetDataN

    subroutine Load( iFileUnit, nCol )
    real(8),allocatable :: Data(:)
    integer iFileUnit, nCol
    integer I,K
    integer byteLen

    integer HeadNo,LenNode
    integer No
    real(8) :: recordTime
    real(8) :: Hertz(3)
    real(8) :: Mrot(3)
    real(8) :: Mtwist(3)
    logical :: is_touching
    logical :: is_slipping
    logical :: is_rolling
    logical :: is_twisting

    type(Nodelink),pointer :: Temp
    type(Nodelink),pointer :: TempH

    allocate(Data(nCol))    
    read(iFileUnit,*) (Data(K),K=1,nCol)
    byteLen = 14

    HeadNo = INT(Data(1))
    LenNode = INT(Data(2))
    !// Load Head Linklist
    Head(HeadNo)%No = LenNode

    Temp => Head(HeadNo)
    if (nCol .GT. 2) then
        do I = 1,(nCol-2)/byteLen
            No = Data(byteLen*(I-1)+3)
            !recordTime = Data(byteLen*(I-1)+4)
            do K = 1,3
                Hertz(K) = Data(byteLen*(I-1)+3+K)
                Mrot(K) = Data(byteLen*(I-1)+6+K)
                Mtwist(K) = Data(byteLen*(I-1)+9+K)
            end do
            is_touching = Data(byteLen*(I-1)+13)
            is_slipping = Data(byteLen*(I-1)+14)
            is_rolling  = Data(byteLen*(I-1)+15)
            is_twisting = Data(byteLen*(I-1)+16)
            !// Load Head Linklist
            allocate(TempH)
            TempH = Nodelink(No,Hertz,Mrot,Mtwist,is_touching,is_slipping,is_rolling,is_twisting,NULL(),NULL())           
            Temp%next => TempH
            TempH%prev => Temp

            Temp => TempH
        end do
    end if       

    deallocate(Data)    
    end subroutine

    subroutine Export( nRow )
    Implicit None
    integer I,J
    integer nRow
    type(Nodelink),pointer :: Temp

    open(13,FILE='export.txt')
    do  J=1,nRow
        write(13,'(2I8)',advance='no') J,Head(J)%No
        if (Head(J)%No .GT. 0) then
            Temp => Head(J)
            do I = 1,Head(J)%No
                write(13,'(I8,2X,9F18.10,2X,4L8)',advance='no') Temp%next%No,Temp%next%Hertz(1),Temp%next%Hertz(2),Temp%next%Hertz(3),Temp%next%Mrot(1),Temp%next%Mrot(2),Temp%next%Mrot(3),Temp%next%Mtwist(1),Temp%next%Mtwist(2),Temp%next%Mtwist(3),Temp%next%is_touching,Temp%next%is_slipping,Temp%next%is_rolling,Temp%next%is_twisting
                Temp => Temp%next
            end do
        end if
        write(13,*)
    end do
    close(13)   

    end subroutine

    Integer Function GetFileN( iFileUnit )
    !// return the number of rows of the file
    Implicit None
    Integer , Intent( IN ) :: iFileUnit
    character( Len = 1 ) :: cDummy
    integer :: ierr
    GetFileN = 0
    Rewind( iFileUnit )
    Do
        Read( iFileUnit , * , ioStat = ierr ) cDummy
        If( ierr /= 0 ) Exit
        GetFileN = GetFileN + 1
    End Do
    Rewind( iFileUnit )
    End Function GetFileN 

    End Module loadFile