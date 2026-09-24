module numbers

implicit none

integer :: Nx, Nt
real(kind=8) :: t, dt, tf, dx, xmin, xmax
real(kind=8), allocatable, dimension(:) :: rhs, u, u_p, x

end module

! -------------------------------------------------------------

program PDE

use numbers
implicit none

integer :: i,n

Nx = 200
tf = 2.0d0
xmin = -1.0d0
xmax = 1.0d0 

allocate(x(0:Nx), rhs(0:Nx), u(0:Nx), u_p(0:Nx))

dx = (xmax - xmin) / dble(Nx)
dt = 0.25d0 * dx

do i = 0, Nx
	x(i) = xmin + dble(i) * dx
end do

Nt = int(tf/dt)

open(1, file='01.dat')

! Condición inicial
t = 0.0d0
u = exp( -x**2 / 0.01d0 )

! Guardar estado inicial
do i = 0, Nx
	write(1,*) x(i), u(i)
end do
write(1,*)

! Integración temporal
do n = 1, Nt

	t = t + dt
	u_p = u
	
	call calcrhs(u_p)
	u = u_p + rhs * dt
		
	call calcrhs(u)
	u = 0.5d0 * ( u_p + u + rhs * dt )

	if ( mod(n,10) == 0 ) then
		write(1,*)		
		write(1,*)
		do i = 0, Nx
			write(1,*) x(i), u(i)
		end do
		print *, n, t
	end if
	
end do

close(1)

end program

! -------------------------------------------------------------

subroutine calcrhs(my_u)

use numbers
implicit none

real(kind=8), dimension(0:Nx), intent(in) :: my_u
integer :: i

do i = 1, Nx-1
	rhs(i) = -0.5d0 * ( my_u(i+1) - my_u(i-1) ) / dx
end do

rhs(0)  = 0.0d0
rhs(Nx) = 0.0d0

end subroutine
