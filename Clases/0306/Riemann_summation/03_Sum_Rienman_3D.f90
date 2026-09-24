program volumen_esfera

implicit none

real(kind=8), allocatable :: x(:,:),y(:,:),f(:,:)
integer :: i,j,Nx,Ny
real(kind=8) :: xmin,xmax,ymin,ymax,dx,dy
real(kind=8) :: R,Integral,Volumen,pi,r2

Nx = 200
Ny = 200

R = 1.0d0

xmin = -R
xmax =  R
ymin = -R
ymax =  R

allocate(x(0:Nx,0:Ny),y(0:Nx,0:Ny),f(0:Nx,0:Ny))

dx = (xmax-xmin)/dble(Nx)
dy = (ymax-ymin)/dble(Ny)

pi = acos(-1.0d0)

! Construcción de la malla
do i=0,Nx
  do j=0,Ny
    x(i,j) = xmin + dble(i)*dx
    y(i,j) = ymin + dble(j)*dy
  end do
end do

! z(x,y) = sqrt(R^2 - x^2 - y^2)
do i=0,Nx
  do j=0,Ny

    r2 = R**2 - x(i,j)**2 - y(i,j)**2

    if (r2 > 0.0d0) then
       f(i,j) = sqrt(r2)
    else
       f(i,j) = 0.0d0
    end if

  end do
end do

! Integración numérica (regla del trapecio en 2D)
Integral = 0.0d0

do i=0,Nx-1
  do j=0,Ny-1
    Integral = Integral + 0.25d0 * &
       ( f(i,j) + f(i+1,j) + f(i,j+1) + f(i+1,j+1) ) * dx * dy
  end do
end do

Volumen = 2.0d0 * Integral

print *, 'Numerical volume =', Volumen
print *, 'Exact volume     =', (4.0d0/3.0d0)*pi*R**3

! Archivo para graficar
open(1,file='Sum_Riemann3D.dat')
do i=0,Nx
  do j=0,Ny
    write(1,*) x(i,j),y(i,j),f(i,j)
  end do
  write(1,*)
end do
close(1)

end program
