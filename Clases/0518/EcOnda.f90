module numbers
  implicit none

  integer :: Nx, Nt
  real(kind=8) :: t, dt, tf, dx, xmin, xmax

  real(kind=8), allocatable :: x(:)

  real(kind=8), allocatable :: phi(:), phi_p(:)
  real(kind=8), allocatable :: psi(:), psi_p(:)
  real(kind=8), allocatable :: pi(:),  pi_p(:)

contains

! ==========================================================
subroutine calcrhs(my_phi,my_psi,my_pi,rhs_phi,rhs_psi,rhs_pi)

  implicit none

  real(kind=8), intent(in)  :: my_phi(0:Nx)
  real(kind=8), intent(in)  :: my_psi(0:Nx)
  real(kind=8), intent(in)  :: my_pi(0:Nx)

  real(kind=8), intent(out) :: rhs_phi(0:Nx)
  real(kind=8), intent(out) :: rhs_psi(0:Nx)
  real(kind=8), intent(out) :: rhs_pi(0:Nx)

  integer :: i

  rhs_phi = my_pi

  ! frontera izquierda
  rhs_psi(0) = 0.5d0 * ( -3.d0*my_pi(0)  &
                       +  4.d0*my_pi(1)  &
                       -       my_pi(2) ) / dx

  rhs_pi(0)  = 0.5d0 * ( -3.d0*my_psi(0) &
                       +  4.d0*my_psi(1) &
                       -       my_psi(2) ) / dx

  ! interior
  do i = 1, Nx-1

    rhs_psi(i) = 0.5d0 * ( my_pi(i+1)  &
                          -my_pi(i-1) ) / dx

    rhs_pi(i)  = 0.5d0 * ( my_psi(i+1) &
                          -my_psi(i-1) ) / dx

  end do

  ! frontera derecha
  rhs_psi(Nx) = 0.5d0 * (  3.d0*my_pi(Nx)   &
                          -4.d0*my_pi(Nx-1) &
                          +      my_pi(Nx-2) ) / dx

  rhs_pi(Nx)  = 0.5d0 * (  3.d0*my_psi(Nx)   &
                          -4.d0*my_psi(Nx-1) &
                          +      my_psi(Nx-2) ) / dx

end subroutine
! ==========================================================

end module

! ==========================================================
program rk2
  use numbers
  implicit none

  integer :: i,n

  real(kind=8), allocatable :: k1_phi(:),k1_psi(:),k1_pi(:)
  real(kind=8), allocatable :: k2_phi(:),k2_psi(:),k2_pi(:)

  xmin = -1.d0
  xmax =  1.d0
  tf   =  2.d0

  Nx = 400

  allocate(x(0:Nx))

  allocate(phi(0:Nx),phi_p(0:Nx))
  allocate(psi(0:Nx),psi_p(0:Nx))
  allocate(pi(0:Nx), pi_p(0:Nx))

  allocate(k1_phi(0:Nx),k1_psi(0:Nx),k1_pi(0:Nx))
  allocate(k2_phi(0:Nx),k2_psi(0:Nx),k2_pi(0:Nx))

  dx = (xmax-xmin)/dble(Nx)
  dt = 0.25d0*dx

  do i=0,Nx
    x(i)=xmin+dble(i)*dx
  end do

  Nt = int(tf/dt)

  ! ======================================================
  ! condiciones iniciales
  ! ======================================================

  phi = exp(-x**2/0.01d0)

  psi = -(2.d0*x/0.01d0)*phi

  pi = 0.d0

  ! ======================================================

  open(1,file='wave1D.dat')

  do i=0,Nx
    write(1,*) x(i),phi(i)
  end do

  ! ======================================================
  ! RK2
  ! ======================================================

  do n=1,Nt

    t = t + dt

    phi_p = phi
    psi_p = psi
    pi_p  = pi

    ! ---------------------------------
    ! k1 = rhs(u^n)
    ! ---------------------------------

    call calcrhs(phi_p,psi_p,pi_p, &
                 k1_phi,k1_psi,k1_pi)

    ! predictor

    phi = phi_p + dt*k1_phi
    psi = psi_p + dt*k1_psi
    pi  = pi_p  + dt*k1_pi

    ! BC
    psi(0) = -0.5d0*(phi(2)-4.d0*phi(1)+3.d0*phi(0))/dx
    pi(0)  =  psi(0)

    psi(Nx)=  0.5d0*(phi(Nx-2)-4.d0*phi(Nx-1)+3.d0*phi(Nx))/dx
    pi(Nx) = -psi(Nx)

    ! ---------------------------------
    ! k2 = rhs(u*)
    ! ---------------------------------

    call calcrhs(phi,psi,pi, &
                 k2_phi,k2_psi,k2_pi)

    ! corrector

    phi = phi_p + 0.5d0*dt*(k1_phi+k2_phi)
    psi = psi_p + 0.5d0*dt*(k1_psi+k2_psi)
    pi  = pi_p  + 0.5d0*dt*(k1_pi +k2_pi )

    ! BC finales
    psi(0) = -0.5d0*(phi(2)-4.d0*phi(1)+3.d0*phi(0))/dx
    pi(0)  =  psi(0)

    psi(Nx)=  0.5d0*(phi(Nx-2)-4.d0*phi(Nx-1)+3.d0*phi(Nx))/dx
    pi(Nx) = -psi(Nx)

    ! salida
    if(mod(n,10)==0) then

      write(1,*)
      write(1,*)

      do i=0,Nx
        write(1,*) x(i),phi(i)
      end do

    end if

  end do

  close(1)

end program
