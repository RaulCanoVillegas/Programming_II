program sum_riemann_2d

implicit none

real(kind=8), allocatable, dimension(:,:) :: x,y,f

integer :: i,j,Nx,Ny

real(kind=8) :: xmin,xmax,dx
real(kind=8) :: ymin,ymax,dy
real(kind=8) :: Integral
real(kind=8) :: pi
real(kind=8) :: r

! Initial conditions
xmin = -5.0d0
xmax =  5.0d0
ymin = -5.0d0
ymax =  5.0d0
Nx = 200
Ny = 200
pi = acos(-1.0d0)

! Allocate memory
allocate(x(0:Nx,0:Ny), y(0:Nx,0:Ny) ,f(0:Nx,0:Ny))	

! Discretization
dx = (xmax - xmin)/dble(Nx)
dy = (ymax - ymin)/dble(Ny)

!Built mesh
do i = 0,Nx
    do j = 0,Ny

        x(i,j) = xmin + dble(i)*dx
        y(i,j) = ymin + dble(j)*dy

    end do
end do

! Circular step function
do i = 0,Nx
    do j = 0,Ny
        r = sqrt(x(i,j)**2 + y(i,j)**2)
        if (r <= 1.0d0) then
            f(i,j) = 1.0d0
        else
            f(i,j) = 0.0d0
        end if
    end do
end do

! Numerical_integrate

Integral = 0.0d0

do i = 0,Nx-1
    do j = 0,Ny-1
        Integral = Integral + &
        0.25d0 * ( f(i,j)     + f(i+1,j) + &
                   f(i,j+1)   + f(i+1,j+1) ) * dx * dy
    end do
end do

print*
print*, 'Numerical integral = ', Integral
print*, 'Exact integral     = ', pi
print*, 'Error              = ', abs(Integral - pi)
print*

! Save
open(1,file='Sum_Riemann2D.dat',status='replace')

do i = 0,Nx
    do j = 0,Ny

        write(1,*) x(i,j),y(i,j),f(i,j)

    end do

    write(1,*)

end do

close(1)

deallocate(x,y,f) 	! memory free

end program
