module numbers

integer :: N,NE
real(kind=8) :: t,dt,tf
real(kind=8) :: s0,i0,r0
real(kind=8) :: beta,gamma,mu
real(kind=8), allocatable, dimension(:) :: k1,k2,k3,k4,rhs,u,u_p

end module

! -------------------------------------------------------------

program rk4

use numbers
implicit none

integer :: i, k
real(kind=8) :: lambda

NE = 3

! Parámetros fijos
gamma = 1.0d0/3.0d0
mu    = 1.0d0/60.0d0

! Condiciones iniciales
s0 = 0.99d0
i0 = 0.01d0
r0 = 0.0d0

tf = 1000.0d0
N  = 10000

allocate(k1(NE), k2(NE), k3(NE), k4(NE), rhs(NE), u(NE), u_p(NE))

dt = tf / dble(N)

open(1, file='2.dat')
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

    write(1,*) k, lambda, t, u(1), u(2), u(3)

  end do
close(1)

end program

! -------------------------------------------------------------

subroutine calcrhs(my_t,my_u)

use numbers
implicit none

real(kind=8), intent(in) :: my_t
real(kind=8), dimension(NE), intent(in) :: my_u

! SIR con demografía
rhs(1) = mu - mu*my_u(1) - beta*my_u(1)*my_u(2)
rhs(2) = beta*my_u(1)*my_u(2) - gamma*my_u(2) - mu*my_u(2)
rhs(3) = gamma*my_u(2) - mu*my_u(3)

end subroutine
