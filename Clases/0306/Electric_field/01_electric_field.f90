program visual_electric_field
implicit none

integer :: i,j,Nx,Ny
real(kind=8) :: xmin,xmax,ymin,ymax,dx,dy
real(kind=8) :: x(0:20,0:20),y(0:20,0:20)
real(kind=8) :: vx(0:20,0:20),vy(0:20,0:20)
real(kind=8) :: r,alpha,escala
character(len=30) :: names

! Parameters
Nx = 20
Ny = 20
xmin = -1.0d0
xmax =  1.0d0
ymin = -1.0d0
ymax =  1.0d0
alpha = 0.2d0
escala = 0.1d0  ! display factor 

! Separation between points
dx = (xmax-xmin)/Nx
dy = (ymax-ymin)/Ny

! Build mesh and calculate field
do i = 0,Nx
    do j = 0,Ny
        x(i,j) = xmin + i*dx
        y(i,j) = ymin + j*dy

        r = sqrt(x(i,j)**2 + y(i,j)**2)		! radial distance

        if ( r > 0.3d0 ) then
            ! Electric radial field normalized for display
            vx(i,j) = escala * x(i,j) / (r**2)  ! decreases as you move away
            vy(i,j) = escala * y(i,j) / (r**2)
        else
            vx(i,j) = 0.0d0
            vy(i,j) = 0.0d0
        end if
    end do
end do

! Guardar en archivo
names = "ElectricField2D.dat"
open(10,file=names	,status="replace")
do i = 0,Nx
    do j = 0,Ny
        write(10,'(4F12.6)') x(i,j), y(i,j), vx(i,j), vy(i,j)
    end do
end do
close(10)

end program 
