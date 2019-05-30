    !********************************************************************
    !     DEMBody 6.0
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
    character(20) :: FileNameB
    character(20) :: FileNameF
    character(25) :: FileNameW
    integer(4) :: result
    integer,external :: systemqq
    real(8) :: temp_bondedWallMeshPoint(3)
    integer :: days,hours,minutes,seconds

    !################         Projectile          ################### 
    !  Update next output time
    Tnext = Tnext + Deltat

    !  Output the position, velocity and radius of the projectile
    open(10,FILE='../Data/Projectile.txt')
    if (Time.LE.Tcrit) then
        write(10,'(20(1XE20.8E4))') Time,(X(K,PP),K=1,3),(Xdot(K,PP),K=1,3),(F(K,PP),K=1,3),(W(K,PP),K=1,3),(FM(K,PP),K=1,3),Energy(PP),(Heat(K,PP),K=1,3)
    end if
    
    if (isBondedWall) then
        !  Output the position, velocity and angular velocity of the Bonded Walls
        open(11,FILE='../Data/BondedWalls.txt')
        if (Time.LE.Tcrit) then
            write(11,'(20(1XE20.8E4))') Time,(bondedWallX(K),K=1,3),(bondedWallXdot(K),K=1,3),(bondedWallW(K),K=1,3),(bondedWallQ(K),K=1,4),(bondedWallF(K),K=1,3),(bondedWallFM(K),K=1,3)
        end if
    end if
    
    if (isBondedTriMeshWall) then
        !  Output the position, velocity and angular velocity of the Bonded TriMesh Walls
        open(12,FILE='../Data/BondedTriMeshWalls.txt')
        if (Time.LE.Tcrit) then
            write(12,'(7(1XE20.8E4))') Time,(bondedTriMeshWallF(K),K=1,3),(bondedTriMeshWallFM(K),K=1,3)
        end if
    end if    
    
    if (isGravBody) then
        !  Output the position, velocity and angular velocity of the Gravity Body
        open(13,FILE='../Data/GravBody.txt')
        if (Time.LE.Tcrit) then
            write(13,'(20(1XE20.8E4))') Time,(gravBodyX(K),K=1,3),(gravBodyXdot(K),K=1,3),(gravBodyW(K),K=1,3),(gravBodyF(K),K=1,3),(gravBodyFM(K),K=1,3),(gravBodyQ(K),K=1,4)
        end if
    end if
    
    if (isSphereBody) then
        !  Output the position, velocity and angular velocity of the Sphere Body
        open(14,FILE='../Data/SphereBody.txt')
        if (Time.LE.Tcrit) then
            do J = 1,sphereBodyNum
                write(14,'(20(1XE20.8E4))') Time,(sphereBodyX(K,J),K=1,3),(sphereBodyXdot(K,J),K=1,3),(sphereBodyW(K,J),K=1,3),(sphereBodyF(K,J),K=1,3),(sphereBodyFM(K,J),K=1,3),(sphereBodyQ(K,J),K=1,4)
            end do
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
            write(Step+1000,*) 'X',',','Y',',','Z',',','U',',','V',',','W',',','R',',','Time',',','EN',',','W:1',',','W:2',',','W:3',',','Slip',',','Roll',',','Twist',',','Q1',',','Q2',',','Q3',',','Q4',',','Cell'
            do  J=1,N
                write(Step+1000,21)X(1,J),',',X(2,J),',',X(3,J),',',Xdot(1,J),',',Xdot(2,J),',',Xdot(3,J),',',R(J),',',Time,',',Energy(J),',',W(1,J),',',W(2,J),',',W(3,J),',',Heat(1,J),',',Heat(2,J),',',Heat(3,J),',',Quaternion(1,J),',',Quaternion(2,J),',',Quaternion(3,J),',',Quaternion(4,J),',',bondTag(J)
            end do
            if (isGravBody) then
                write(Step+1000,21)gravBodyX(1),',',gravBodyX(2),',',gravBodyX(3),',',gravBodyXdot(1),',',gravBodyXdot(2),',',gravBodyXdot(3),',',gravBodyR,',',Time,',',0.0,',',gravBodyW(1),',',gravBodyW(2),',',gravBodyW(3),',',0.0,',',0.0,',',0.0,',',gravBodyQ(1),',',gravBodyQ(2),',',gravBodyQ(3),',',gravBodyQ(4),',',0
            end if
            if (isSphereBody) then
                do J=1,sphereBodyNum
                    write(Step+1000,21)sphereBodyX(1,J),',',sphereBodyX(2,J),',',sphereBodyX(3,J),',',sphereBodyXdot(1,J),',',sphereBodyXdot(2,J),',',sphereBodyXdot(3,J),',',sphereBodyR(J),',',Time,',',0.0,',',sphereBodyW(1,J),',',sphereBodyW(2,J),',',sphereBodyW(3,J),',',0.0,',',0.0,',',0.0,',',sphereBodyQ(1,J),',',sphereBodyQ(2,J),',',sphereBodyQ(3,J),',',sphereBodyQ(4,J),',',0
                end do
            end if            
21          format (E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,I8)
            close(Step+1000)
        else
            open(Step+1000,FILE=FileNameX)
            write(Step+1000,*) 'X',',','Y',',','Z',',','U',',','V',',','W',',','R',',','Time',',','EN',',','W:1',',','W:2',',','W:3',',','Slip',',','Roll',',','Twist',',','Cell'
            do  J=1,N
                write(Step+1000,22)X(1,J),',',X(2,J),',',X(3,J),',',Xdot(1,J),',',Xdot(2,J),',',Xdot(3,J),',',R(J),',',Time,',',Energy(J),',',W(1,J),',',W(2,J),',',W(3,J),',',Heat(1,J),',',Heat(2,J),',',Heat(3,J),',',bondTag(J)
            end do
            if (isGravBody) then
                write(Step+1000,22)gravBodyX(1),',',gravBodyX(2),',',gravBodyX(3),',',gravBodyXdot(1),',',gravBodyXdot(2),',',gravBodyXdot(3),',',gravBodyR,',',Time,',',Energy(J),',',gravBodyW(1),',',gravBodyW(2),',',gravBodyW(3),',',X0(1,J),',',X0(2,J),',',X0(3,J),',',0
            end if
22          format (E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,I8)
            close(Step+1000)
        end if
        
        !  Output the position, velocity and force of all assemblies
        write(FileNameB,'(I4)') Step+1000
        FileNameB = '../Data/'//FileNameB
        FileNameB = trim(FileNameB)//'B.csv'
        open(Step+1000,FILE=FileNameB)
        
        write(Step+1000,*) 'X',',','Y',',','Z',',','U',',','V',',','W',',','Time',',','W:1',',','W:2',',','W:3',',','Q1',',','Q2',',','Q3',',','Q4'
        do  J=1,bondN
            write(Step+1000,23)bondX(1,J),',',bondX(2,J),',',bondX(3,J),',',bondXdot(1,J),',',bondXdot(2,J),',',bondXdot(3,J),',',Time,',',bondW(1,J),',',bondW(2,J),',',bondW(3,J),',',bondQ(1,J),',',bondQ(2,J),',',bondQ(3,J),',',bondQ(4,J)
        end do
23      format (E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6,A2,E15.6)
        close(Step+1000)

#ifdef historyOutput        
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
                    write(Step+1000,'(I8,2X,9(2XF18.8),2X,4L8)',advance='no') Temp%next%No,Temp%next%Hertz(1),Temp%next%Hertz(2),Temp%next%Hertz(3),Temp%next%Mrot(1),Temp%next%Mrot(2),Temp%next%Mrot(3),Temp%next%Mtwist(1),Temp%next%Mtwist(2),Temp%next%Mtwist(3),Temp%next%is_touching,Temp%next%is_slipping,Temp%next%is_rolling,Temp%next%is_twisting
                    Temp => Temp%next
                end do
            end if
            write(Step+1000,*)
        end do
        close(Step+1000)
#endif        
        
        !  Output the Bonded Wall Meshfile
        if (isBondedWall) then
            write(FileNameW,'(I4)') Step+1000
#ifdef Linux
            FileNameW = '../Data/Wall/'//FileNameW
            FileNameW = trim(FileNameW)//'W.vtk'            
            result = systemqq("cp ../Input/bondedWallMesh.vtk "//FileNameW)
#elif Windows
            FileNameW = '..\Data\Wall\'//FileNameW
            FileNameW = trim(FileNameW)//'W.vtk'   
            result = systemqq("copy ..\Input\bondedWallMesh.vtk "//FileNameW)
#endif
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

#ifdef openmp
    T2 = omp_get_wtime()
#endif
    days = int((T2-T1)/3600/24)
    hours = int(((T2-T1)-days*3600*24)/3600)
    minutes = int(((T2-T1)-days*3600*24-hours*3600)/60)
    seconds = int((T2-T1)-days*3600*24-hours*3600-minutes*60)
    
    write(*,"(A11,I4,1X,A4,I4,1X,A4,I4,1X,A4,I4,1X,A4)") "time cost: ", days, "days", hours, "hours", minutes, "minutes", seconds, "seconds"        
    write(*,*) "refresh frequency: ", dble(refreshNum)/(Tcrit/Dt+1)

    stop
    end