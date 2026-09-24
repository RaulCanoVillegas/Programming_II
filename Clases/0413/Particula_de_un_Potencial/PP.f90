module numbers

implicit none

integer :: N, NE
real(kind=8) :: t, dt, tf
real(kind=8) :: G, M, a, e
real(kind=8) :: x0, y0, z0, vx0, vy0, vz0
real(kind=8), allocatable :: k1(:), k2(:), k3(:), k4(:)
real(kind=8), allocatable :: rhs(:), u(:), u_p(:)

end module

! -------------------------------------------------------------

program Particula_en_un_Potencial

use numbers
implicit none

integer :: i

NE = 6

! Constantes físicas
G  = 6.6743d-11
M  = 1.98847d30
a  = 5.906423d12
e  = 0.2488

! Condiciones iniciales (perihelio)
x0  = a*(1.0d0 - e**2)/(1.0d0 + e)
y0  = 0.0d0
z0  = 0.0d0

vx0 = 0.0d0
vy0 = sqrt(G*M/(a*(1.0d0 - e**2))) * (1.0d0 - e**2)
vz0 = 0.0d0

! Tiempo
tf = 60.0d0*60.0d0*24.0d0*365.0d0*300.0d0   ! 300 años
dt = 100.0d0
N  = int(tf/dt)

allocate(k1(NE), k2(NE), k3(NE), k4(NE), rhs(NE), u(NE), u_p(NE))

open(1, file='Pluto.dat')

! Inicialización
t = 0.0d0
u(1) = x0
u(2) = y0
u(3) = z0
u(4) = vx0
u(5) = vy0
u(6) = vz0

! Integración RK4
do i = 1, N

  u_p = u

  call calcrhs(t, u_p)
  k1 = rhs

  call calcrhs(t + 0.5d0*dt, u_p + 0.5d0*dt*k1)
  k2 = rhs

  call calcrhs(t + 0.5d0*dt, u_p + 0.5d0*dt*k2)
  k3 = rhs

  call calcrhs(t + dt, u_p + dt*k3)
  k4 = rhs

  u = u_p + (dt/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)

  t = t + dt

  write(1,*) t, u(1), u(2), u(3)

end do

close(1)

end program

! -------------------------------------------------------------

subroutine calcrhs(my_t, my_u)

use numbers
implicit none

real(kind=8), intent(in) :: my_t
real(kind=8), dimension(NE), intent(in) :: my_u

real(kind=8) :: x, y, z, r

! Posiciones actuales
x = my_u(1)
y = my_u(2)
z = my_u(3)

r = sqrt(x*x + y*y + z*z)

! Ecuaciones
rhs(1) = my_u(4)
rhs(2) = my_u(5)
rhs(3) = my_u(6)

rhs(4) = -G*M*x/(r**(3/2))
rhs(5) = -G*M*y/(r**(3/2))
rhs(6) = -G*M*z/(r**(3/2))

end subroutine
