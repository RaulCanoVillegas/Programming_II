
module numbers

integer :: Nx,Nt
real(kind=8) :: t,dt,tf,dx,xmin,xmax
real(kind=8), allocatable, dimension(:) :: x
real(kind=8), allocatable, dimension(:) :: rhs_phi, phi, phi_p
real(kind=8), allocatable, dimension(:) :: rhs_pi, pi, pi_p
real(kind=8), allocatable, dimension(:) :: rhs_psi, psi, psi_p
real(kind=8) :: a,c

end module

! -------------------

program Wave

use numbers
implicit none

integer :: i,n

tf = 2.0d0
Nx  = 400
xmin = -1.0d0
xmax = 1.0d0
c = 0.25d0

allocate(x(0:Nx), rhs_phi(0:Nx), phi(0:Nx), phi_p(0:Nx))
allocate(rhs_pi(0:Nx), pi(0:Nx), pi_p(0:Nx))
allocate(rhs_psi(0:Nx), psi(0:Nx), psi_p(0:Nx))

dx = ( xmax-xmin) / dble(Nx)
dt = c*dx

do i=0, Nx
  x(i) = xmin + dble(i) * dx !Dominio discreto
end do

Nt = int(tf/dt)

! Initial conditions
t = 0.0d0


open(1,file='Wave')

phi = exp(-(x**2.0d0)/0.010d0)
psi = phi * (-(2d0*x)/0.01d0)
pi = 0.0d0

do i=0, Nx
  write(1,*)  x(i), phi(i)
end do



do n=1,Nt

  t    = t + dt
  phi_p  = phi
  psi_p  = psi
  pi_p  = pi


  ! ----------------
  call calcrhs( phi_p, psi_p, pi_p  )
    phi  = phi_p+ rhs_phi*dt
    psi  = psi_p + rhs_psi*dt
    pi  = pi_p + rhs_pi*dt

  psi(0) = 0.5d0 * (phi(2) - 4d0*phi(1) + 3* phi(0) ) / dx
  pi(0) = psi(0)
  psi(Nx) = 0.5d0 * (phi(Nx-2) - 4d0*phi(Nx-1) + 3* phi(Nx) ) / dx
  pi(Nx) = psi(Nx)

  call calcrhs( phi, psi, pi  )
    phi  = 0.5d0*(phi_p + phi + rhs_phi*dt)
    psi  = 0.5d0*(psi_p + psi + rhs_psi*dt)
    pi  = 0.5d0*(pi_p + pi + rhs_pi*dt)

  psi(0) = 0.5d0 * (phi(2) - 4d0*phi(1) + 3* phi(0) ) / dx
  pi(0) = psi(0)
  psi(Nx) = 0.5d0 * (phi(Nx-2) - 4d0*phi(Nx-1) + 3* phi(Nx) ) / dx
  pi(Nx) = psi(Nx)
  ! ----------------
  if( mod(n,10).eq.0 ) then
    write(1,*)
    write(1,*)
    do i=0, Nx
      write(1,*)  x(i), phi(i)
    end do
  end if
end do
close(1)

end program

! ----------------------
subroutine calcrhs(my_phi, my_psi, my_pi )

use numbers
implicit none


real(kind=8), dimension(0:Nx), intent(in) :: my_psi, my_phi, my_pi
integer :: i

do i=1, Nx-1
  rhs_phi(i) = my_pi(i)
  rhs_psi(i) = 0.5d0 * (my_pi(i+1) - my_pi(i-1)) / dx 
  rhs_pi(i) = 0.5d0 * (my_psi(i+1) - my_psi(i-1)) / dx 
end do

rhs_phi(0) = my_pi(0)
rhs_psi(0) = 0.5d0 * (my_pi(1) - my_pi(Nx-1)) / dx 
rhs_pi(0) = 0.5d0 * (my_psi(1) - my_psi(Nx-1)) / dx 

rhs_phi(Nx) = my_pi(0)
rhs_psi(Nx) = 0.5d0 * (my_pi(1) - my_pi(Nx-1)) / dx 
rhs_pi(Nx) = 0.5d0 * (my_psi(1) - my_psi(Nx-1)) / dx 


end subroutine
