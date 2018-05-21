    program main

    integer :: I,J,K,L
    real(8) :: RadiusBound  !mesh grid
    real(8) :: R,width  !particle params
    integer :: number,iterate
    real(8) :: r_rand,theta_rand,psi_rand
    real(8) :: tol
    real(8) :: gamma
    real(8),parameter :: PI = 3.1415926 

    integer :: Tag
    integer :: Flag
    integer :: nRow
    real(8) :: point(3)
    real(8),allocatable :: points(:,:)
    real(8) :: dx,dist(3)
    real(8) :: radius
    real(8) :: velocity(3)

    !--- Params initialize
    RadiusBound = 150.0

    R = 4.5
    width = 0.5
    number = 20000
    tol = 0.1
    iterate = 100
    gamma = 0
    shooting = 15
    allocate(points(4,number))
    
    points = 0.0    

    open(10,FILE='particles.txt')
    call random_seed()
    do I = 1,number
        write(*,*) I
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
                    
                    if (abs(point(1)).GE.D .OR. abs(point(2)).GE.D) then
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

    end

