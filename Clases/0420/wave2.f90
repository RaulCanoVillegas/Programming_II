module numbers
  integer :: Nx, Nt
  real(kind=8) :: t, dt, tf, dx, xmin, xmax
  real(kind=8), allocatable, dimension(:) :: x
  real(kind=8), allocatable, dimension(:) :: rhs_phi, phi, phi_p
  real(kind=8), allocatable, dimension(:) :: rhs_pi,  pi,  pi_p
  real(kind=8), allocatable, dimension(:) :: rhs_psi, psi, psi_p
end module

! -------------------
program rk2
  use numbers
  implicit none
  integer :: i, n

  xmin = -1.0d0
  xmax =  1.0d0
  tf   =  2.0d0
  Nx   =  400

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

  open(1, file='wave.dat')

  ! Condiciones iniciales
  t   = 0.0d0
  phi = exp(-x**2 / 0.01d0)
  psi = -(2.0d0 * x * phi) / 0.01d0
  pi  = 0.0d0

  do i = 0, Nx
    write(1,*) x(i), t, phi(i)
  end do

  ! Bucle temporal
  do n = 1, Nt

    t     = t + dt
    phi_p = phi
    psi_p = psi
    pi_p  = pi

    ! --- Paso 1 del RK2 ---
    call calcrhs(phi_p, psi_p, pi_p)

    phi = phi_p + rhs_phi * dt
    psi = psi_p + rhs_psi * dt
    pi  = pi_p  + rhs_pi  * dt

    ! Condiciones de frontera Dirichlet
    phi(0)  = 0.0d0
    phi(Nx) = 0.0d0
    pi(0)   = 0.0d0
    pi(Nx)  = 0.0d0
    psi(0)  = (-3*phi(0) + 4*phi(1)  - phi(2))   / (2*dx)
    psi(Nx) = ( 3*phi(Nx)- 4*phi(Nx-1)+ phi(Nx-2))/ (2*dx)

    ! --- Paso 2 del RK2 ---
    call calcrhs(phi, psi, pi)

    phi = 0.5d0*(phi_p + phi + rhs_phi * dt)
    psi = 0.5d0*(psi_p + psi + rhs_psi * dt)
    pi  = 0.5d0*(pi_p  + pi  + rhs_pi  * dt)

    ! Condiciones de frontera Dirichlet
    phi(0)  = 0.0d0
    phi(Nx) = 0.0d0
    pi(0)   = 0.0d0
    pi(Nx)  = 0.0d0
    psi(0)  = (-3*phi(0)  + 4*phi(1)   - phi(2))    / (2*dx)
    psi(Nx) = ( 3*phi(Nx) - 4*phi(Nx-1)+ phi(Nx-2)) / (2*dx)

    ! Escribir cada 10 pasos
    if (mod(n,10) == 0) then
      write(1,*)
      write(1,*)
      do i = 0, Nx
        write(1,*) x(i), t, phi(i)
      end do
    end if

  end do

  close(1)
end program

! ----------------------
subroutine calcrhs(my_phi, my_psi, my_pi)
  use numbers
  implicit none
  real(kind=8), dimension(0:Nx), intent(in) :: my_phi, my_psi, my_pi
  integer :: i

  do i = 1, Nx-1
    rhs_psi(i) = 0.5d0*(my_pi(i+1)  - my_pi(i-1))  / dx
    rhs_pi(i)  = 0.5d0*(my_psi(i+1) - my_psi(i-1)) / dx
  end do

  ! Fronteras periódicas para rhs (interior)
  rhs_psi(0)  = 0.5d0*(my_pi(1)  - my_pi(Nx-1))  / dx
  rhs_psi(Nx) = rhs_psi(0)
  rhs_pi(0)   = 0.5d0*(my_psi(1) - my_psi(Nx-1)) / dx
  rhs_pi(Nx)  = rhs_pi(0)

  rhs_phi = my_pi

end subroutine
