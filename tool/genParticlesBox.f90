    !********************************************************************
    !     genParticles
    !     ************
    !
    !     DEM & NBODY CODE
    !     ----------------
    !
    !     Developed by Bin Cheng, Tsinghua.
    !     .................................
    !
    !********************************************************************
    program main
    implicit none
    
    integer :: I,J,K,L
    real(8) :: D(2),H(2),T(2)  !mesh grid
    real(8) :: R,width  !particle params
    integer :: number,iterate,shooting
    real(8) :: tol
    real(8) :: gamma
    real(8) :: mass,surface
    real(8) :: Rho

    integer :: Tag
    integer :: Flag
    integer :: nRow
    real(8) :: point(3)
    real(8) :: dx,dist(3)
    real(8) :: radius
    real(8) :: velocity(3)
    real(8),allocatable :: points(:,:)

    integer :: r_law_n
    integer,allocatable :: r_law(:)
    real(8) :: r_law_step
    real(8) :: temp(4)
    integer :: step
    
    integer :: Distribution
    real(8),external :: power_law
    real(8) :: power
    
    Distribution = 0
    !################         Input         ###################
    !--- Mesh Region
    D = (/-20,20/)
    H = (/-20,20/)
    T = (/0,60/)
    
    !--- Particle Params
    R = 1.5
    width = 0.5
    Rho = 2.5
    
    !--- Particle Number
    number = 2000
    tol = 0.01
    iterate = 200
    gamma = 0.0
    shooting = 100
    allocate(points(4,number))
    
    !--- Power Law
    Distribution = 0
    power = 2.8
    r_law_n = 100
    r_law_step = width*2.0/REAL(r_law_n)
    allocate(r_law(r_law_n))
    r_law = 0
    
    !################         Part 1          ###################
    !--- Generation
    points = 0.0    
    
    open(10,FILE='particles.txt')
    call random_seed()
    do I = 1,number
        if (MODULO(I,1000)==0) then
            write(*,*) I
        end if
        Tag = 1
        do while (Tag < iterate)
            
            if (Distribution == 1) then
                call random_number(radius)
                radius = power_law(R-width, R+width, radius, power) 
            else
                call random_number(radius)
                radius = radius*width*2.0-width+R
            end if
            
            call random_number(point(1))
            point(1) = (D(2)-D(1))*point(1)+D(1)
            call random_number(point(2))
            point(2) = (H(2)-H(1))*point(2)+H(1)
            call random_number(point(3))
            point(3) = (T(2)-T(1))*point(3)+T(1)
    
            Flag = 1
            if (I > 1) then
                do J = 1,I-1 !All points
                    do K = 1,3
                        dist(K) = point(K) - points(K,J)
                    end do
                    dx = sqrt(dist(1)**2 + dist(2)**2 + dist(3)**2)
                    dx = dx - radius - points(4,J)
                    if (dx < tol) then
                        Flag = 0
                        exit
                    end if
                end do
            end if
            
            !---shooting method (only five steps)
            if (Flag == 0) then
                do L = 1,shooting
                    do K = 1,3
                        point(K) = point(K) + dist(K)*dx*2.0
                    end do
                    
                    if (point(1).LE.D(1) .OR. point(1).GE.D(2) .OR. point(2).LE.H(1) .OR. point(2).GE.H(2) .OR. point(3).LE.T(1) .OR. point(3).GE.T(2)) then
                        exit
                    end if
    
                    Flag = 1
                    if (I > 1) then
                        do J = 1,I-1 !All points
                            do K = 1,3
                                dist(K) = point(K) - points(K,J)
                            end do
                            dx = sqrt(dist(1)**2 + dist(2)**2 + dist(3)**2)
                            dx = dx - radius - points(4,J)
                            if (dx < tol) then
                                Flag = 0
                                exit
                            end if
                        end do
                    end if
                    
                    if (Flag == 1) then
                        exit
                    end if
                end do
            end if
    
            if (Flag == 1) then
                call random_number(velocity(1))
                velocity(1) = point(2)*gamma + (velocity(1)*0.002-0.001)
                call random_number(velocity(2))
                velocity(2) = (velocity(2)*0.002-0.001)
                call random_number(velocity(3))
                velocity(3) = (velocity(3)*0.002-0.001)
                points(1,I) = point(1)
                points(2,I) = point(2)
                points(3,I) = point(3)
                points(4,I) = radius
                write(10,'(7F15.5)') (points(K,I),K=1,4),(velocity(K),K=1,3)
                exit
            end if
    
            Tag = Tag + 1
        end do
    
    end do
    close(10)
    
    !################         Part 2          ###################
    mass = 0.0D0
    surface = 0.0D0
    open(10,FILE='particles.txt')
    do I = 1,number
        read(10,*) (temp(K),K=1,4)
        step = FLOOR((temp(4)-(R-width))/r_law_step)+1
        if (step .LE. r_law_n) then
            r_law(step) = r_law(step) + 1
        end if
        mass = mass + 4./3.*3.1415926*temp(4)**3
        surface = surface + 3.1415926*temp(4)**2
    end do
    close(10)
    
    open(10,FILE='r_law.txt')
    do I = 1,r_law_n
        write(10,*) r_law(I)
    end do
    close(10)
    write(*,*) mass,surface
    
    !################         Part 3          ###################
    open(10,FILE='particles.txt')
    open(11,FILE='output.txt')
    open(12,FILE='paraview.csv')
    write(12,*) 'X',',','Y',',','Z',',','U',',','V',',','W',',','R' 
    do I = 1,number
        read (10,*) (point(K),K=1,3),radius,(velocity(K),K=1,3)
        write(11,'(12F15.5)') 4./3.*3.1415926*radius**3*Rho,0.4*4./3.*3.1415926*radius**5*Rho,(point(K),K=1,3),(velocity(K),K=1,3),0.0,0.0,0.0,radius
        write(12,'(F12.5,A2,F12.5,A2,F12.5,A2,F12.5,A2,F12.5,A2,F12.5,A2,F12.5)') point(1),',',point(2),',',point(3),',',velocity(1),',',velocity(2),',',velocity(3),',',radius
    end do
    end
    
    function power_law(r_min, r_max, y, gamma)
    implicit none
    real(8) :: r_min,r_max,y,gamma
    reaL(8) :: power_law
    
    power_law = ((r_max**(-gamma+1) - r_min**(-gamma+1))*y  + r_min**(-gamma+1.0))**(1.0/(-gamma + 1.0))
    
    return
end

