    !********************************************************************
    !     DEMBody 3.0
    !     ***********
    !
    !     Output and data save.
    !     ---------------------
    !
    !********************************************************************
    subroutine output()
    
    use global
    use omp_lib
    implicit none
    
    integer :: I,J,K
    type(Nodelink),pointer :: Temp
    character(20) :: FileNameX
    character(20) :: FileNameF
    character(25) :: FileNameW
    integer(4) :: result
    integer,external :: systemqq
    real(8) :: temp_bondedWallMeshPoint(3)

    !################         Projectile          ################### 
    !  Update next output time
    Tnext = Tnext + Deltat

    !  Output the position, velocity and radius of the projectile
    open(10,FILE='../Data/Projectile.txt')
    if (Time.LE.Tcrit) then
        write(10,'(14E20.8)') Time,(X(K,N),K=1,3),(Xdot(K,N),K=1,3),(F(K,N),K=1,3),(W(K,N),K=1,3),Energy(N)
    end if
    
    if (isBondedWall) then
        !  Output the position, velocity and angular velocity of the Bonded Walls
        open(11,FILE='../Data/BondedWalls.txt')
        if (Time.LE.Tcrit) then
            write(11,'(20E15.7)') Time,(bondedWallX(K),K=1,3),(bondedWallXdot(K),K=1,3),(bondedWallW(K),K=1,3),(bondedWallQ(K),K=1,4),(bondedWallF(K),K=1,3),(bondedWallFM(K),K=1,3)
        end if
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
            write(Step+1000,*) 'X',',','Y',',','Z',',','U',',','V',',','W',',','F:1',',','F:2',',','F:3',',','R',',','Time',',','EN',',','W:1',',','W:2',',','W:3',',','X0',',','Y0',',','Z0',',','color'
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
        
        !  Output the Bonded Wall Meshfile
        if (isBondedWall) then
            write(FileNameW,'(I4)') Step+1000
            FileNameW = '../Data/Wall/'//FileNameW
            FileNameW = trim(FileNameW)//'W.vtk'
            
            result = systemqq("cp ../Input/bondedWallMesh.vtk "//FileNameW)
            open(Step+1000,FILE=FileNameW,position='Append')
            do J=1,size(bondedWallMeshPoint,2)
                do K = 1,3
                    temp_bondedWallMeshPoint(K) = bondedWallX(K) + bondedWallMatB(K,1)*bondedWallMeshPoint(1,J) + bondedWallMatB(K,2)*bondedWallMeshPoint(2,J) + bondedWallMatB(K,3)*bondedWallMeshPoint(3,J)
                end do
                write(Step+1000,'(3F10.4)') (temp_bondedWallMeshPoint(K),K=1,3)
            end do  
            close(Step+1000)
        end if
        
        Step = Step + 1          
    end if

    !  Check termination time.
    if (Time<Tcrit) return

    T2 = omp_get_wtime()
    write(*,*) "time cost: ", (T2-T1)

    stop
    end