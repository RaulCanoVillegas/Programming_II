module numbers

integer :: N,NE
real(kind=8) :: t,dt,tf,x0,y0
real(kind=8) :: a,b,c,d
real(kind=8), allocatable, dimension(:) :: k1,k2,k3,k4,rhs,u,u_p

end module

! -------------------------------------------------------------

program rk4

use numbers
implicit none

integer :: i

NE = 2

! Parámetros del modelo depredador-presa
a = 1.25d0
b = 0.7d0
c = 0.4d0
d = 0.25d0

! Condiciones iniciales
x0 = 3.0d0
y0 = 2.0d0

! Intervalo de tiempo
tf = 100.0d0
N  = 1000

allocate(k1(1:NE),k2(1:NE),k3(1:NE),k4(1:NE),rhs(1:NE),u(1:NE),u_p(1:NE))

dt = tf / dble(N)

! Inicialización
t = 0.0d0
u(1) = x0
u(2) = y0

open(1,file='1.dat')
write(1,*) t, u(1), u(2)

do i=1,N

  u_p = u

! --------------------------- RK4 ---------------------------
  call calcrhs(t, u_p)
  k1 = rhs

  call calcrhs(t + 0.5d0*dt, u_p + 0.5d0*dt*k1)
  k2 = rhs

  call calcrhs(t + 0.5d0*dt, u_p + 0.5d0*dt*k2)
  k3 = rhs

  call calcrhs(t + dt, u_p + dt*k3)
  k4 = rhs

  u = u_p + (dt/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
! -------------------------------------------------------------

  t = t + dt

  write(1,*) t, u(1), u(2)

end do

close(1)

end program

! -------------------------------------------------------------

subroutine calcrhs(my_t,my_u)

use numbers
implicit none

real(kind=8), intent(in) :: my_t
real(kind=8), dimension(NE), intent(in) :: my_u

! Sistema Lotka-Volterra (introducir ecuaciones)
rhs(1) = a*my_u(1) - b*my_u(1)*my_u(2)
rhs(2) = -c*my_u(2) + d*my_u(1)*my_u(2)

end subroutine
