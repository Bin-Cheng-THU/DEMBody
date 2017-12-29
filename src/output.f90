    !********************************************************************
    !     DEMBody 2.0
    !     ***********
    !
    !     Output and data save.
    !     ---------------------
    !
    !********************************************************************
    subroutine output()

    use global
    implicit none
    
    integer :: I,J,K
    type(Nodelink),pointer :: Temp
    character(20) :: FileNameX
    character(20) :: FileNameF

    !################         Projectile          ################### 
    !  Update next output time
    Tnext = Tnext + Deltat

    !  Output the position, velocity and radius of the projectile
    open(10,FILE='../Data/Projectile.txt')
    if (Time.LE.Tcrit) then
        write(10,'(14E20.8)') Time,(X(K,N),K=1,3),(Xdot(K,N),K=1,3),(F(K,N),K=1,3),(W(K,N),K=1,3),Energy(N)
    end if

    if (Time .GE. CheckPointTnext) then
        !  Update next check point time
        CheckPointTnext = CheckPointTnext + CheckPointDt

        !  Output the position, velocity and force of all particles
        write(FileNameX,'(I4)') Step+1000
        FileNameX = '../Data/'//FileNameX
        FileNameX = trim(FileNameX)//'X.csv'

        if (isQuaternion) then
            open(Step+1000,FILE=FileNameX)
            write(Step+1000,*) 'X',',','Y',',','Z',',','U',',','V',',','W',',','F:1',',','F:2',',','F:3',',','R',',','Time',',','EN',',','W:1',',','W:2',',','W:3',',','X0',',','Y0',',','Z0',',','Q1',',','Q2',',','Q3',',','Q4'
            do  J=1,N
                write(Step+1000,21)X(1,J),',',X(2,J),',',X(3,J),',',Xdot(1,J),',',Xdot(2,J),',',Xdot(3,J),',',F(1,J),',',F(2,J),',',F(3,J),',',R(J),',',Time,',',Energy(J),',',W(1,J),',',W(2,J),',',W(3,J),',',X0(1,J),',',X0(2,J),',',X0(3,J),',',Quaternion(1,J),',',Quaternion(2,J),',',Quaternion(3,J),',',Quaternion(4,J)
            end do
21          format (F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15)      
            close(Step+1000)
        else
            open(Step+1000,FILE=FileNameX)
            write(Step+1000,*) 'X',',','Y',',','Z',',','U',',','V',',','W',',','F:1',',','F:2',',','F:3',',','R',',','Time',',','EN',',','W:1',',','W:2',',','W:3',',','X0',',','Y0',',','Z0'
            do  J=1,N
                write(Step+1000,22)X(1,J),',',X(2,J),',',X(3,J),',',Xdot(1,J),',',Xdot(2,J),',',Xdot(3,J),',',F(1,J),',',F(2,J),',',F(3,J),',',R(J),',',Time,',',Energy(J),',',W(1,J),',',W(2,J),',',W(3,J),',',X0(1,J),',',X0(2,J),',',X0(3,J)
            end do
22          format (F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15,A2,F30.15)      
            close(Step+1000)
        end if

        !  Output the tangential contact history of all interactions
        write(FileNameF,'(I4)') Step+1000
        FileNameF = '../Data/'//FileNameF
        FileNameF = trim(FileNameF)//'F.txt'

        open(Step+1000,FILE=FileNameF)
        do  J=1,N
            write(Step+1000,'(2I8)',advance='no') J,Head(J)%No
            if (Head(J)%No .GT. 0) then
                Temp => Head(J)
                do I = 1,Head(J)%No
                    write(Step+1000,'(I8,2X,6F15.5,2X,3I8)',advance='no') Temp%next%No,Temp%next%Hertz(1),Temp%next%Hertz(2),Temp%next%Hertz(3),Temp%next%Mrot(1),Temp%next%Mrot(2),Temp%next%Mrot(3),Temp%next%is_touching,Temp%next%is_slipping,Temp%next%is_rolling
                    Temp => Temp%next
                end do
            end if
            write(Step+1000,*)
        end do
        close(Step+1000)

        Step = Step + 1          
    end if

    !  Check termination time.
    if (Time<Tcrit) return

    call CPU_TIME(T2)
    write(*,*) "time cost: ", (T2-T1)/20/3600

    stop
    end