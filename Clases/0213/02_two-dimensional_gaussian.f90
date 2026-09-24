program twodimiensional_gaussian

implicit none

real(kind=8), allocatable, dimension(:,:) :: x,y,f
integer :: i,j,Nx,Ny
real(kind=8) :: xmin, xmax, dx, ymin, ymax, dy

! Initial conditions
xmin = -5.0
xmax = 5.0
ymin = -5.0
ymax = 5.0
Nx = 100
Ny = 100

! Allocate memory
allocate(x(0:Nx,0:Ny),y(0:Nx,0:Ny),f(0:Nx,0:Ny))

! Separation between points
dx = ( xmax - xmin ) / dble(Nx)
dy = ( ymax - ymin ) / dble(Ny)

! Build (x, y)
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

! Function
f = exp(-x**2 -y**2)

open(1,file='two-dimensional.dat')
do i=0,Nx
  do j=0,Ny
    write(1,*) x(i,j),y(i,j),f(i,j)
  end do
  write(1,*)
end do
close(1)

end program
