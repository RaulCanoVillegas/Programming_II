program electric_field_anim

implicit none

integer :: i,j,k,Nx,Ny,Nt
real(kind=8) :: xmin,xmax,ymin,ymax,dx,dy
real(kind=8) :: t,dt,pi
real(kind=8) :: r,escala

real(kind=8) :: x(0:20,0:20),y(0:20,0:20)
real(kind=8) :: vx(0:20,0:20),vy(0:20,0:20)

! Parameters

Nx = 20
Ny = 20
Nt = 100
xmin = -1.0d0
xmax =  1.0d0
ymin = -1.0d0
ymax =  1.0d0
pi = acos(-1.0d0)
escala = 0.1d0

! Separation between points
dx = (xmax-xmin)/Nx
dy = (ymax-ymin)/Ny
dt = 2.0d0*pi/dble(Nt)

! Built mesh
do i = 0,Nx
    do j = 0,Ny

        x(i,j) = xmin + i*dx
        y(i,j) = ymin + j*dy

    end do
end do

! Save
open(10,file='animated_field.dat',status='replace')
do k = 0,Nt
    t = k*dt
    do i = 0,Nx
        do j = 0,Ny
            r = sqrt(x(i,j)**2 + y(i,j)**2)

            if (r > 0.15d0) then
                vx(i,j) = escala * x(i,j)/(r**2) * cos(t)
                vy(i,j) = escala * y(i,j)/(r**2) * cos(t)
            else
                vx(i,j) = 0.0d0
                vy(i,j) = 0.0d0
            end if
            write(10,*) x(i,j),y(i,j),vx(i,j),vy(i,j)
        end do
        write(10,*)
    end do
    write(10,*)
    write(10,*)
end do

close(10)

end program
