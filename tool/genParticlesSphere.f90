    program main
    
    implicit none
    integer :: I,J,K,L
    real(8) :: RadiusBound  !mesh grid
    real(8) :: R,width  !particle params
    integer :: number,iterate
    real(8) :: r_rand,theta_rand,psi_rand
    real(8) :: tol
    real(8) :: gamma
    real(8),parameter :: PI = 3.1415926 
    integer :: shooting
    
    integer :: Tag
    integer :: Flag
    integer :: nRow
    real(8) :: point(3)
    real(8),allocatable :: points(:,:)
    real(8) :: dx,dist(3)
    real(8) :: radius
    real(8) :: velocity(3)
    real(8) :: rho

    !--- Params initialize
    RadiusBound = 5.5

    R = 0.2
    width = 0.05
    rho = 2.5
    number = 3000
    tol = 0.1
    iterate = 200
    gamma = 0
    shooting = 50
    allocate(points(4,number))
    
    points = 0.0    

    open(10,FILE='particles.txt')
    call random_seed()
    do I = 1,number
        if (mod(I,50)==0) then
            write(*,*) I
        end if
        Tag = 1
        do while (Tag < iterate)
            call random_number(radius)
            radius = 2*width*radius-width+R 
            
            call random_number(r_rand)
            r_rand = (RadiusBound-radius)*r_rand
            call random_number(theta_rand)
            theta_rand = PI*theta_rand
            call random_number(psi_rand)
            psi_rand = 2.0*PI*psi_rand
            
            point(1) = r_rand*sin(theta_rand)*cos(psi_rand)
            point(2) = r_rand*sin(theta_rand)*sin(psi_rand)
            point(3) = r_rand*cos(theta_rand)

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

