program discrete_domain

implicit none
real(kind=8), allocatable, dimension(:,:) :: x,y,vx,vy
real(kind=8) :: xmin,xmax,ymin,ymax,dx,dy,dt,t
integer :: i,j,Nx,Ny,k
real(kind=8) pi

! Initial conditions
pi = acos(-1.0d0)
xmin = -1.0d0
xmax =  1.0d0
ymin = -1.0d0
ymax =  1.0d0
Nx   = 20
Ny   = 20

! Allocate memory
allocate(x(0:Nx,0:Ny),y(0:Nx,0:Ny),vx(0:Nx,0:Ny),vy(0:Nx,0:Ny))

! Discrete domain
dx   = (xmax - xmin) / dble(Nx)
dy   = (ymax - ymin) / dble(Ny)
do i=0,Nx
  do j=0,Ny
    x(i,j) = xmin + dble(i) * dx
    y(i,j) = ymin + dble(j) * dy
  end do
end do

! Time
t  = 0.0d0
dt = 4.0d0*pi / 100.0d0

! Save
open(1,file='vector_field.dat')

DO k=0,100

  vx = - 0.2 * ( y / sqrt( x**2 + y**2 ) ) * cos(t)	! component x
  vy =   0.2 * ( x / sqrt( x**2 + y**2 ) ) * cos(t)	! component y

  ! Output
  write(1,*) '# Time =',t
  do i=0,Nx
    do j=0,Ny
      write(1,*) x(i,j),y(i,j),vx(i,j),vy(i,j)
    end do
    write(1,*)
  end do
  write(1,*)
  write(1,*)
  t = t + dt

END DO

close(1)

deallocate(x,y,vx,vy)

end program

