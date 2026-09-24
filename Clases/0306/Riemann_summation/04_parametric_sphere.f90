program esfera_parametrica

implicit none

integer :: i,j,Nt,Np
real(kind=8) :: theta,phi,dtheta,dphi
real(kind=8) :: x,y,z,R,pi

Nt = 100
Np = 100

R = 1.0d0
pi = acos(-1.0d0)

dtheta = pi/dble(Nt)
dphi   = 2.0d0*pi/dble(Np)

open(1,file='sphere.dat')

do i=0,Nt
  theta = dble(i)*dtheta
  do j=0,Np

    phi = dble(j)*dphi

    x = R*sin(theta)*cos(phi)
    y = R*sin(theta)*sin(phi)
    z = R*cos(theta)

    write(1,*) x,y,z

  end do
  write(1,*)
end do

close(1)

end program esfera_parametrica
