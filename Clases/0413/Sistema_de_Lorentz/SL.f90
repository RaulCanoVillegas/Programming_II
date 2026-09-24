module numbers

implicit none

integer :: N, NE
real(kind=8) :: t, dt, tf
real(kind=8) :: x0, y0, z0
real(kind=8) :: a, b, c
real(kind=8), allocatable, dimension(:) :: k1, k2, k3, k4, rhs, u, u_p

end module

! -------------------------------------------------------------

program rk4

use numbers
implicit none

integer :: i

NE = 3

! Parámetros del sistema de Lorenz
a = 15.0d0
b = 30.0d0
c = 3.0d0

! Condiciones iniciales
x0 = 7.0d0
y0 = 6.0d0
z0 = 5.0d0

! Tiempo
tf = 100.0d0
N  = 100000        ! dt = 0.001

allocate(k1(NE), k2(NE), k3(NE), k4(NE), rhs(NE), u(NE), u_p(NE))

dt = tf / dble(N)

open(1, file='2.dat')

! Inicialización
t = 0.0d0
u(1) = x0
u(2) = y0
u(3) = z0

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

! Sistema de Lorenz
rhs(1) = a * (my_u(2) - my_u(1))
rhs(2) = b * my_u(1) - my_u(2) - my_u(1) * my_u(3)
rhs(3) = my_u(1) * my_u(2) - c * my_u(3)

end subroutine
