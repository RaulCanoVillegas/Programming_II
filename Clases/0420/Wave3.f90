module numbers
  implicit none
  integer :: Nx, Nt
  real(kind=8) :: t, dt, tf, dx, xmin, xmax
  real(kind=8), allocatable :: x(:)
  real(kind=8), allocatable :: rhs_phi(:), phi(:), phi_p(:)
  real(kind=8), allocatable :: rhs_pi(:),  pi(:),  pi_p(:)
  real(kind=8), allocatable :: rhs_psi(:), psi(:), psi_p(:)
end module

! -------------------
program rk2
  use numbers
  implicit none
  integer :: i, n

  ! Parámetros
  Nx   = 200
  xmin = -1.0d0
  xmax =  1.0d0
  tf   =  0.25d0

  ! Memoria
  allocate(x(0:Nx))
  allocate(rhs_phi(0:Nx), phi(0:Nx), phi_p(0:Nx))
  allocate(rhs_pi(0:Nx),  pi(0:Nx),  pi_p(0:Nx))
  allocate(rhs_psi(0:Nx), psi(0:Nx), psi_p(0:Nx))

  dx = (xmax - xmin) / dble(Nx)
  dt = 0.25d0 * dx

  do i = 0, Nx
    x(i) = xmin + dble(i)*dx
  end do

  Nt = int(tf / dt)

  open(1, file='wave_200.dat')

  ! Condiciones iniciales (gaussiana)
  t   = 0.0d0
  phi = exp(-x**2 / 0.01d0)
  psi = -(2.0d0 * x * phi) / 0.01d0
  pi  = 0.0d0

  do i = 0, Nx
    write(1,*) x(i), t, phi(i)
  end do

  ! Evolución temporal
  do n = 1, Nt

    t     = t + dt
    phi_p = phi
    psi_p = psi
    pi_p  = pi

    ! ---- Paso 1 ----
    call calcrhs(phi_p, psi_p, pi_p)

    phi = phi_p + rhs_phi * dt
    psi = psi_p + rhs_psi * dt
    pi  = pi_p  + rhs_pi  * dt

    call boundary_conditions

    ! ---- Paso 2 ----
    call calcrhs(phi, psi, pi)

    phi = 0.5d0*(phi_p + phi + rhs_phi * dt)
    psi = 0.5d0*(psi_p + psi + rhs_psi * dt)
    pi  = 0.5d0*(pi_p  + pi  + rhs_pi  * dt)

    call boundary_conditions

    if (mod(n,10) == 0) then
      write(1,*)
      write(1,*)
      do i = 0, Nx
        write(1,*) x(i), t, phi(i)
      end do
    end if

  end do

  close(1)
  write(*,*) 'Terminado Nx=200'

end program

! ----------------------
subroutine calcrhs(my_phi, my_psi, my_pi)
  use numbers
  implicit none
  real(kind=8), intent(in) :: my_phi(0:Nx), my_psi(0:Nx), my_pi(0:Nx)
  integer :: i

  ! Diferencias centradas (NO periódicas)
  do i = 1, Nx-1
    rhs_psi(i) = (my_pi(i+1)  - my_pi(i-1))  / (2.0d0*dx)
    rhs_pi(i)  = (my_psi(i+1) - my_psi(i-1)) / (2.0d0*dx)
  end do

  rhs_phi = my_pi

end subroutine

! ----------------------
subroutine boundary_conditions
  use numbers
  implicit none

  ! Fronteras tipo cuerda fija + salida
  phi(0)  = 0.0d0
  phi(Nx) = 0.0d0

  pi(0)   = psi(0)
  pi(Nx)  = -psi(Nx)

  psi(0)  = (-3*phi(0)  + 4*phi(1)   - phi(2))    / (2.0d0*dx)
  psi(Nx) = ( 3*phi(Nx) - 4*phi(Nx-1)+ phi(Nx-2)) / (2.0d0*dx)

end subroutine
