module numbers
implicit none

integer :: Nx, Nt
real(kind=8) :: t, dt, tf, dx, xmin, xmax

real(kind=8), allocatable :: x(:)

real(kind=8), allocatable :: phi(:), psi(:), pi(:)
real(kind=8), allocatable :: phi_p(:), psi_p(:), pi_p(:)

real(kind=8), allocatable :: k1_phi(:), k2_phi(:), k3_phi(:), k4_phi(:)
real(kind=8), allocatable :: k1_psi(:), k2_psi(:), k3_psi(:), k4_psi(:)
real(kind=8), allocatable :: k1_pi(:),  k2_pi(:),  k3_pi(:),  k4_pi(:)

real(kind=8), allocatable :: rhs_phi(:), rhs_psi(:), rhs_pi(:)

end module

! -------------------------------------------------------------

program wave
use numbers
implicit none

integer :: i, n

! Dominio espacial
Nx = 200
xmin = -1.0d0
xmax = 1.0d0
dx = (xmax - xmin) / dble(Nx-1)

! Tiempo
tf = 2.0d0
Nt = 2000
dt = tf / dble(Nt)

! Memoria
allocate(x(Nx))

allocate(phi(Nx), psi(Nx), pi(Nx))
allocate(phi_p(Nx), psi_p(Nx), pi_p(Nx))

allocate(rhs_phi(Nx), rhs_psi(Nx), rhs_pi(Nx))

allocate(k1_phi(Nx), k2_phi(Nx), k3_phi(Nx), k4_phi(Nx))
allocate(k1_psi(Nx), k2_psi(Nx), k3_psi(Nx), k4_psi(Nx))
allocate(k1_pi(Nx),  k2_pi(Nx),  k3_pi(Nx),  k4_pi(Nx))

! Mallado
do i = 1, Nx
    x(i) = xmin + (i-1)*dx
end do

! -------------------------------------------------------------
! 🔥 CONDICIONES INICIALES: dos gaussianas en las fronteras
! -------------------------------------------------------------

do i = 1, Nx
    phi(i) = exp(-100.0d0*(x(i)+0.8d0)**2) + exp(-100.0d0*(x(i)-0.8d0)**2)
end do

! Velocidad inicial → ambas viajan hacia el centro
do i = 2, Nx-1
    pi(i) = - (phi(i+1) - phi(i-1)) / (2.0d0*dx)
end do

pi(1) = 0.0d0
pi(Nx) = 0.0d0

! psi consistente
do i = 2, Nx-1
    psi(i) = (phi(i+1) - phi(i-1)) / (2.0d0*dx)
end do

psi(1) = 0.0d0
psi(Nx) = 0.0d0

t = 0.0d0

open(1,file='wave.dat')

! -------------------------------------------------------------
! Evolución temporal (RK4)
! -------------------------------------------------------------
do n = 1, Nt

    phi_p = phi
    psi_p = psi
    pi_p  = pi

    ! k1
    call calcrhs(phi_p, psi_p, pi_p)
    k1_phi = rhs_phi
    k1_psi = rhs_psi
    k1_pi  = rhs_pi

    ! k2
    call calcrhs(phi_p + 0.5d0*dt*k1_phi, psi_p + 0.5d0*dt*k1_psi, pi_p + 0.5d0*dt*k1_pi)
    k2_phi = rhs_phi
    k2_psi = rhs_psi
    k2_pi  = rhs_pi

    ! k3
    call calcrhs(phi_p + 0.5d0*dt*k2_phi, psi_p + 0.5d0*dt*k2_psi, pi_p + 0.5d0*dt*k2_pi)
    k3_phi = rhs_phi
    k3_psi = rhs_psi
    k3_pi  = rhs_pi

    ! k4
    call calcrhs(phi_p + dt*k3_phi, psi_p + dt*k3_psi, pi_p + dt*k3_pi)
    k4_phi = rhs_phi
    k4_psi = rhs_psi
    k4_pi  = rhs_pi

    ! Actualización
    phi = phi_p + (dt/6.0d0)*(k1_phi + 2.0d0*k2_phi + 2.0d0*k3_phi + k4_phi)
    psi = psi_p + (dt/6.0d0)*(k1_psi + 2.0d0*k2_psi + 2.0d0*k3_psi + k4_psi)
    pi  = pi_p  + (dt/6.0d0)*(k1_pi  + 2.0d0*k2_pi  + 2.0d0*k3_pi  + k4_pi)

    ! ---------------------------------------------------------
    ! 🔴 Fronteras tipo cuerda fija
    ! ---------------------------------------------------------
    phi(1)  = 0.0d0
    phi(Nx) = 0.0d0

    psi(1)  = 0.0d0
    psi(Nx) = 0.0d0

    pi(1)  = 0.0d0
    pi(Nx) = 0.0d0

    t = t + dt

    ! Guardar
    if (mod(n,20) == 0) then
        do i = 1, Nx
            write(1,*) x(i), phi(i)
        end do
        write(1,*)
        write(1,*)
    end if

end do

close(1)

end program

! -------------------------------------------------------------

subroutine calcrhs(phi_in, psi_in, pi_in)

use numbers
implicit none

real(kind=8), intent(in) :: phi_in(Nx), psi_in(Nx), pi_in(Nx)
integer :: i

do i = 2, Nx-1

    rhs_phi(i) = pi_in(i)

    rhs_psi(i) = (pi_in(i+1) - pi_in(i-1)) / (2.0d0*dx)

    rhs_pi(i)  = (psi_in(i+1) - psi_in(i-1)) / (2.0d0*dx)

end do

rhs_phi(1)  = 0.0d0
rhs_phi(Nx) = 0.0d0

rhs_psi(1)  = 0.0d0
rhs_psi(Nx) = 0.0d0

rhs_pi(1)   = 0.0d0
rhs_pi(Nx)  = 0.0d0

end subroutine
