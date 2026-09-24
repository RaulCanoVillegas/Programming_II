program animated_gradient

implicit none

integer :: i,j,n,Nx,Ny,Nt
real(8) :: xmin,xmax,ymin,ymax,dx,dy,dt
real(8) :: x,y,t
real(8) :: dfdx, dfdy
real(8), parameter :: pi = 4.d0*atan(1.d0)

! Domain
xmin = -pi
xmax =  pi
ymin = -pi
ymax =  pi

! Variables
Nx = 30
Ny = 30
Nt = 60

! Separation between points
dx = (xmax - xmin)/(Nx-1)
dy = (ymax - ymin)/(Ny-1)
dt = 2*pi/Nt			! frames

! Save
open(10,file='gradient.dat')

! Temporal cycle
do n = 0, Nt-1 
    t = n*dt
    do i = 0, Nx-1
        x = xmin + i*dx
        do j = 0, Ny-1
            y = ymin + j*dy
            dfdx =  3.d0*cos(3.d0*x)*cos(4.d0*y)*cos(t)		! partial derivative with respect to x
            dfdy = -4.d0*sin(3.d0*x)*sin(4.d0*y)*cos(t)		! partial derivative with respect to y	
            write(10,*) x,y,dfdx,dfdy
        end do
        write(10,*)
    end do
    write(10,*)
    write(10,*)
end do
close(10)

end program 
