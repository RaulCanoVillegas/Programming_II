module numbers
  implicit none
  integer :: Nx, Ny, Nt
  real(kind=8) :: t, dt, tf, dx, dy, xmin, xmax, ymin, ymax
  real(kind=8), allocatable, dimension(:)   :: x, y
  real(kind=8), allocatable, dimension(:,:) :: phi,   phi_p
  real(kind=8), allocatable, dimension(:,:) :: psi_x, psi_x_p
  real(kind=8), allocatable, dimension(:,:) :: psi_y, psi_y_p
  real(kind=8), allocatable, dimension(:,:) :: pi_f,  pi_p
  real(kind=8), allocatable, dimension(:,:) :: rhs_phi
  real(kind=8), allocatable, dimension(:,:) :: rhs_psi_x
  real(kind=8), allocatable, dimension(:,:) :: rhs_psi_y
  real(kind=8), allocatable, dimension(:,:) :: rhs_pi
end module numbers

! -------------------------------------------------------
program water_drop
  use numbers
  implicit none
  integer :: i, j, n
  real(kind=8) :: sigma2

  ! --- domain setup ---
  xmin = -2.0d0;  xmax = 2.0d0
  ymin = -2.0d0;  ymax = 2.0d0
  tf   = 4.0d0
  Nx   = 200;  Ny = 200

  allocate( x(0:Nx), y(0:Ny) )
  allocate( phi(0:Nx,0:Ny),   phi_p(0:Nx,0:Ny)   )
  allocate( psi_x(0:Nx,0:Ny), psi_x_p(0:Nx,0:Ny) )
  allocate( psi_y(0:Nx,0:Ny), psi_y_p(0:Nx,0:Ny) )
  allocate( pi_f(0:Nx,0:Ny),  pi_p(0:Nx,0:Ny)    )
  allocate( rhs_phi(0:Nx,0:Ny)   )
  allocate( rhs_psi_x(0:Nx,0:Ny) )
  allocate( rhs_psi_y(0:Nx,0:Ny) )
  allocate( rhs_pi(0:Nx,0:Ny)    )

  dx = (xmax - xmin) / dble(Nx)
  dy = (ymax - ymin) / dble(Ny)
  dt = 0.25d0 * min(dx, dy)       ! CFL condition

  do i = 0, Nx
    x(i) = xmin + dble(i)*dx
  end do
  do j = 0, Ny
    y(j) = ymin + dble(j)*dy
  end do

  Nt = int(tf / dt)

  ! --- initial condition: Gaussian "drop" centered at origin ---
  sigma2 = 0.1d0
  do j = 0, Ny
    do i = 0, Nx
      phi(i,j)   =  exp( -(x(i)**2 + y(j)**2) / sigma2 )
      psi_x(i,j) = -(2.0d0*x(i)/sigma2) * phi(i,j)
      psi_y(i,j) = -(2.0d0*y(j)/sigma2) * phi(i,j)
      pi_f(i,j)  =  0.0d0
    end do
  end do

  ! --- output: write initial data ---
  open(1, file='water_drop.dat')
  t = 0.0d0
  call write_slice(1, t)

  ! =====================  RK2 time loop  ======================
  do n = 1, Nt

    t = t + dt

    ! --- save current level ---
    phi_p   = phi
    psi_x_p = psi_x
    psi_y_p = psi_y
    pi_p    = pi_f

    ! ---- Stage 1: u* = u^n + dt * RHS(u^n) ----
    call calcrhs(phi_p, psi_x_p, psi_y_p, pi_p)
    phi   = phi_p   + rhs_phi   * dt
    psi_x = psi_x_p + rhs_psi_x * dt
    psi_y = psi_y_p + rhs_psi_y * dt
    pi_f  = pi_p    + rhs_pi    * dt
    call apply_bc(phi, psi_x, psi_y, pi_f)

    ! ---- Stage 2: u^{n+1} = 0.5*(u^n + u* + dt*RHS(u*)) ----
    call calcrhs(phi, psi_x, psi_y, pi_f)
    phi   = 0.5d0*(phi_p   + phi   + rhs_phi   * dt)
    psi_x = 0.5d0*(psi_x_p + psi_x + rhs_psi_x * dt)
    psi_y = 0.5d0*(psi_y_p + psi_y + rhs_psi_y * dt)
    pi_f  = 0.5d0*(pi_p    + pi_f  + rhs_pi    * dt)
    call apply_bc(phi, psi_x, psi_y, pi_f)

    ! --- output every 20 steps (slice along y = Ny/2) ---
    if (mod(n, 20) == 0) then
      call write_slice(1, t)
    end if

  end do

  close(1)
  write(*,*) 'Done.  Output written to water_drop.dat'
end program water_drop

! -------------------------------------------------------
!  RHS of the first-order system:
!    phi_t  = pi
!    psi_x_t = pi_x   (= d/dx pi)
!    psi_y_t = pi_y   (= d/dy pi)
!    pi_t   = psi_x_x + psi_y_y  (= div(grad phi))
! -------------------------------------------------------
subroutine calcrhs(my_phi, my_psi_x, my_psi_y, my_pi)
  use numbers
  implicit none
  real(kind=8), intent(in), dimension(0:Nx,0:Ny) :: my_phi, my_psi_x, my_psi_y, my_pi
  integer :: i, j

  do j = 1, Ny-1
    do i = 1, Nx-1
      rhs_phi(i,j)   = my_pi(i,j)
      rhs_psi_x(i,j) = 0.5d0*(my_pi(i+1,j) - my_pi(i-1,j)) / dx
      rhs_psi_y(i,j) = 0.5d0*(my_pi(i,j+1) - my_pi(i,j-1)) / dy
      rhs_pi(i,j)    = 0.5d0*(my_psi_x(i+1,j) - my_psi_x(i-1,j)) / dx &
                     + 0.5d0*(my_psi_y(i,j+1) - my_psi_y(i,j-1)) / dy
    end do
  end do

  ! boundaries set to zero (handled in apply_bc)
  rhs_phi(0,:)    = 0.0d0;  rhs_phi(Nx,:)   = 0.0d0
  rhs_phi(:,0)    = 0.0d0;  rhs_phi(:,Ny)   = 0.0d0
  rhs_psi_x(0,:)  = 0.0d0;  rhs_psi_x(Nx,:) = 0.0d0
  rhs_psi_x(:,0)  = 0.0d0;  rhs_psi_x(:,Ny) = 0.0d0
  rhs_psi_y(0,:)  = 0.0d0;  rhs_psi_y(Nx,:) = 0.0d0
  rhs_psi_y(:,0)  = 0.0d0;  rhs_psi_y(:,Ny) = 0.0d0
  rhs_pi(0,:)     = 0.0d0;  rhs_pi(Nx,:)    = 0.0d0
  rhs_pi(:,0)     = 0.0d0;  rhs_pi(:,Ny)    = 0.0d0
end subroutine calcrhs

! -------------------------------------------------------
!  Outgoing (Sommerfeld-like) radiation boundary conditions
!  at all four walls.  One-sided 2nd-order differences.
! -------------------------------------------------------
subroutine apply_bc(my_phi, my_psi_x, my_psi_y, my_pi)
  use numbers
  implicit none
  real(kind=8), intent(inout), dimension(0:Nx,0:Ny) :: my_phi, my_psi_x, my_psi_y, my_pi
  integer :: i, j

  ! --- x = xmin wall (i=0) ---
  do j = 0, Ny
    my_psi_x(0,j) = -0.5d0*( my_phi(2,j) - 4.0d0*my_phi(1,j) + 3.0d0*my_phi(0,j) ) / dx
    my_pi(0,j)    =  my_psi_x(0,j)
  end do

  ! --- x = xmax wall (i=Nx) ---
  do j = 0, Ny
    my_psi_x(Nx,j) =  0.5d0*( my_phi(Nx-2,j) - 4.0d0*my_phi(Nx-1,j) + 3.0d0*my_phi(Nx,j) ) / dx
    my_pi(Nx,j)    = -my_psi_x(Nx,j)
  end do

  ! --- y = ymin wall (j=0) ---
  do i = 0, Nx
    my_psi_y(i,0) = -0.5d0*( my_phi(i,2) - 4.0d0*my_phi(i,1) + 3.0d0*my_phi(i,0) ) / dy
    my_pi(i,0)    =  my_psi_y(i,0)
  end do

  ! --- y = ymax wall (j=Ny) ---
  do i = 0, Nx
    my_psi_y(i,Ny) =  0.5d0*( my_phi(i,Ny-2) - 4.0d0*my_phi(i,Ny-1) + 3.0d0*my_phi(i,Ny) ) / dy
    my_pi(i,Ny)    = -my_psi_y(i,Ny)
  end do
end subroutine apply_bc

! -------------------------------------------------------
!  Write a 1D slice (along x at j=Ny/2) to the output file.
!  Format compatible with gnuplot's 'index' blocks.
! -------------------------------------------------------
subroutine write_slice(unit_no, time)
  use numbers
  implicit none
  integer,       intent(in) :: unit_no
  real(kind=8),  intent(in) :: time
  integer :: i, jmid

  jmid = Ny / 2
  write(unit_no, *)
  write(unit_no, *)
  write(unit_no, '(A,F10.5)') '# t = ', time
  do i = 0, Nx
    write(unit_no, *) x(i), phi(i, jmid)
  end do
end subroutine write_slice
