module numbers

integer :: N,NE
real(kind=8) :: t,dt,tf,x0,y0,k,m
real(kind=8) :: F0,b,Omega
real(kind=8) :: xexacta,xerror
real(kind=8), allocatable, dimension(:) :: k1,k2,k3,k4,rhs,u,u_p

end module

! -------------------

program rk4

use numbers
implicit none

integer :: i

NE = 2
k  = 2.0d0
m  = 0.5d0	
x0 = 1.0d0
y0 = 0.0d0

F0    = 0.2d0
b     = 0.05d0
Omega = 0.5d0

tf = 1.0d1
N  = 200

allocate(k1(1:NE),k2(1:NE),k3(1:NE),k4(1:NE),rhs(1:NE),u(1:NE),u_p(1:NE))

dt = tf / dble(N)

! Initial conditions
t = 0.0d0
u(1) = x0
u(2) = y0

open(1,file='rk4_200.dat')
write(1,*) t, u(1),u(2)

do i=1,N

  u_p = u

  ! ----------- RK4 -------------
  call calcrhs(t, u_p)
  k1 = rhs

  call calcrhs(t + 0.5d0*dt, u_p + 0.5d0*dt*k1)
  k2 = rhs

  call calcrhs(t + 0.5d0*dt, u_p + 0.5d0*dt*k2)
  k3 = rhs

  call calcrhs(t + dt, u_p + dt*k3)
  k4 = rhs

  u = u_p + (dt/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
  ! -----------------------------

  t = t + dt

  write(1,*) t, u(1),u(2)

end do

close(1)

end program

! ----------------------

subroutine calcrhs(my_t,my_u)

use numbers
implicit none

real(kind=8), intent(in) :: my_t
real(kind=8), dimension(NE), intent(in) :: my_u

rhs(1) = my_u(2)
rhs(2) = -k/m * my_u(1) - b/m * my_u(2) + F0/m * sin(Omega * my_t)

end subroutine
