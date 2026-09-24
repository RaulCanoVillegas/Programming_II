
module numbers

integer :: Nx, Nt
real(kind=8) :: t, dt, tf, dx, xmin, xmax
real(kind=8), allocatable, dimension(:) :: x
real(kind=8), allocatable, dimension(:) :: phi, phi_p, phi_pp

! phi_old = phi_p
! phi_new = phi_pp

end module

! -------------------

program wave_transparent

use numbers
implicit none

integer :: i, n
real(kind=8) :: sigma, r

xmin = -1.0d0
xmax =  1.0d0
tf   =  3.0d0
Nx   =  400

allocate( x(0:Nx) )
allocate( phi(0:Nx), phi_pp(0:Nx), phi_p(0:Nx) )

dx = (xmax - xmin) / dble(Nx)
dt = 0.5d0 * dx     ! CFL = dt/dx = 0.5

r = dt / dx         ! numero de Courant

do i = 0, Nx
   x(i) = xmin + dble(i)*dx
end do
Nt = int(tf / dt)

open(1, file='onda_transparente.dat')

! ---- Condicion inicial ----
sigma = 0.1d0
t     = 0.0d0
phi     = exp(-x**2 / sigma**2)
phi_p = phi    ! phi en t=-dt => velocidad inicial cero

! Escribir frame inicial
do i = 0, Nx
   write(1,*) x(i), t, phi(i)
end do

! ---- Bucle temporal ----
do n = 1, Nt

   t = t + dt

   ! --- Puntos interiores ---
   do i = 1, Nx-1
      phi_pp(i) = 2.0d0*phi(i) - phi_p(i) &
                   + r**2 * (phi(i+1) - 2.0d0*phi(i) + phi(i-1))
   end do

   ! --- Frontera Izquierda ---
   phi_pp(0) = phi(1) + (r - 1.0d0)/(r + 1.0d0) * (phi_pp(1) - phi(0))

   ! --- Frontera Derecha ---
   phi_pp(Nx) = phi(Nx-1) + (r - 1.0d0)/(r + 1.0d0) * (phi_pp(Nx-1) - phi(Nx))

   ! tiempos
   phi_p = phi
   phi     = phi_pp

   if (mod(n,10) .eq. 0) then
      write(1,*)
      write(1,*)
      do i = 0, Nx
         write(1,*) x(i), t, phi(i)
      end do
   end if

end do

close(1)

end program
