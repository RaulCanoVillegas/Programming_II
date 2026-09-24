program gaussian2D_time

implicit none

real(kind=8), allocatable, dimension(:,:) :: x,y,f
integer :: i,j,Nx,Ny,Nt,k
real(kind=8) :: xmin, xmax, dx, ymin, ymax, dy, t, dt, pi

! Initial conditions
xmin = -5.0
xmax = 5.0
ymin = -5.0
ymax = 5.0
Nx = 100
Ny = 100
Nt = 100
pi = acos(-1.0d0)

! Allocate memory
allocate(x(0:Nx,0:Ny),y(0:Nx,0:Ny),f(0:Nx,0:Ny))

! Separation between points
dx = ( xmax - xmin ) / dble(Nx)
dy = ( ymax - ymin ) / dble(Ny)
dt = 2.0d0 * pi / dble(Nt)		! time data

! Build (variables)
do i=0,Nx
  do j=0,Ny
    x(i,j) = xmin + dble(i) * dx
  end do
end do
do i=0,Nx
  do j=0,Ny
    y(i,j) = ymin + dble(j) * dy
  end do
end do

! Save
t = 0.0
open(1,file='gaussian2D_time.dat')
do k=0,Nt 				! time cycle
  f = exp(-x**2 -y**2) * cos(t)		! time evolution 
  do i=0,Nx
    do j=0,Ny
      write(1,*) x(i,j),y(i,j),f(i,j)
    end do
    write(1,*)
  end do
  write(1,*)
  write(1,*)
  t = t + dt				! advance time
end do
close(1)

end program
