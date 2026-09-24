module numbers
implicit none
integer :: Nx, Nt
real(kind=8) :: t, dt, tf, dx, xmin, xmax, c
real(kind=8), allocatable :: phi(:), psi(:), phi_new(:), psi_new(:), x(:)
end module

! -------------------------------------------------------------
program wave
use numbers
implicit none
integer :: i, n

Nx = 200
tf = 2.0d0
xmin = -1.0d0
xmax = 1.0d0
c = 1.0d0

allocate(x(0:Nx), phi(0:Nx), psi(0:Nx), phi_new(0:Nx), psi_new(0:Nx))

dx = (xmax - xmin) / dble(Nx)
dt = 0.4d0 * dx / c   ! CFL condition

do i = 0, Nx
    x(i) = xmin + dble(i)*dx
end do

Nt = int(tf/dt)

open(1,file='wave.dat')

! -------------------------------
! Condiciones iniciales: DOS gaussianas
! -------------------------------
do i = 0, Nx
    phi(i) = exp(-(x(i)+0.5d0)**2/0.01d0) + exp(-(x(i)-0.5d0)**2/0.01d0)
    psi(i) = 0.0d0
end do

! Guardar inicial
call save()

! -------------------------------
! Evolución temporal
! -------------------------------
do n = 1, Nt

    call calcrhs(phi, psi, phi_new, psi_new)

    phi = phi_new
    psi = psi_new

    if (mod(n,10)==0) then
        call save()
        print*, n, t
    end if

    t = t + dt
end do

close(1)

end program

! -------------------------------------------------------------
subroutine calcrhs(phi, psi, phi_new, psi_new)
use numbers
implicit none
real(kind=8), intent(in)  :: phi(0:Nx), psi(0:Nx)
real(kind=8), intent(out) :: phi_new(0:Nx), psi_new(0:Nx)
integer :: i
real(kind=8) :: d2phi

! Interior
do i = 1, Nx-1
    d2phi = (phi(i+1) - 2.0d0*phi(i) + phi(i-1)) / dx**2

    phi_new(i) = phi(i) + dt * psi(i)
    psi_new(i) = psi(i) + dt * c**2 * d2phi
end do

! -------------------------------
! Condiciones de frontera (rebote)
! -------------------------------
phi_new(0) = 0.0d0
phi_new(Nx) = 0.0d0

psi_new(0) = 0.0d0
psi_new(Nx) = 0.0d0

end subroutine

! -------------------------------------------------------------
subroutine save()
use numbers
implicit none
integer :: i

do i = 0, Nx
    write(1,*) x(i), phi(i)
end do
write(1,*)
write(1,*)

end subroutine
