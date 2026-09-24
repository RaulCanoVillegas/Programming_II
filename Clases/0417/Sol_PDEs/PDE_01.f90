module numbers
  implicit none
  integer  :: Nx, Nt
  real(kind=8) :: t, dt, tf, dx, xmin, xmax, c_wave
  real(kind=8), allocatable, dimension(:) :: rhs_u, rhs_v, u, v, u_p, v_p, x
end module

! -------------------------------------------------------------
program WaveCollision
  use numbers
  implicit none
  integer :: i, n

  Nx      = 400
  tf      = 3.2d0
  xmin    = -1.0d0
  xmax    =  1.0d0
  c_wave  =  1.0d0

  allocate(x(0:Nx), rhs_u(0:Nx), rhs_v(0:Nx), &
           u(0:Nx), v(0:Nx), u_p(0:Nx), v_p(0:Nx))

  dx = (xmax - xmin) / dble(Nx)
  dt = 0.4d0 * dx / c_wave   ! condición CFL

  do i = 0, Nx
    x(i) = xmin + dble(i) * dx
  end do

  Nt = int(tf / dt)

  open(1, file='wave_collision.dat')

  ! ── Condición inicial ──────────────────────────────────────
  ! Dos gaussianas: izquierda viaja hacia +x, derecha hacia -x
  t = 0.0d0
  call set_initial(1.0d0)

  ! Guardar estado inicial
  do i = 0, Nx
    write(1,*) t, x(i), u(i)
  end do
  write(1,*)

  ! ── Integración temporal (RK2 predictor-corrector) ─────────
  do n = 1, Nt
    t   = t + dt
    u_p = u
    v_p = v

    ! Paso predictor (Euler)
    call calcrhs(u_p, v_p)
    u = u_p + rhs_u * dt
    v = v_p + rhs_v * dt

    ! Paso corrector (promedio Heun)
    call calcrhs(u, v)
    u = 0.5d0 * (u_p + u + rhs_u * dt)
    v = 0.5d0 * (v_p + v + rhs_v * dt)

    ! Guardar cada 10 pasos
    if (mod(n, 10) == 0) then
      write(1,*)
      write(1,*)
      do i = 0, Nx
        write(1,*) t, x(i), u(i)
      end do
      print *, 'paso:', n, '  t =', t
    end if

  end do

  close(1)
  print *, 'Listo. Datos en wave_collision.dat'

contains

  subroutine set_initial(amp)
    real(kind=8), intent(in) :: amp
    real(kind=8) :: sig2, x0, xi, gL, gR
    integer :: i
    sig2 = 0.015d0   ! anchura de las gaussianas
    x0   = 0.55d0    ! posición inicial (±x0)
    do i = 0, Nx
      xi = x(i)
      gL = amp * exp(-(xi + x0)**2 / sig2)   ! gaussiana izquierda
      gR = amp * exp(-(xi - x0)**2 / sig2)   ! gaussiana derecha
      u(i) = gL + gR
      ! Velocidad inicial: gL viaja hacia +x, gR hacia -x
      ! v = du/dt = c * d(gL)/dx - c * d(gR)/dx
      v(i) =  c_wave * gL * 2.0d0*(xi + x0) / sig2 &
             -c_wave * gR * 2.0d0*(xi - x0) / sig2
    end do
  end subroutine

end program

! -------------------------------------------------------------
subroutine calcrhs(my_u, my_v)
  use numbers
  implicit none
  real(kind=8), dimension(0:Nx), intent(in) :: my_u, my_v
  integer :: i

  ! du/dt = v
  rhs_u(0:Nx) = my_v(0:Nx)

  ! dv/dt = c^2 * d²u/dx²  (diferencias centradas de 2do orden)
  do i = 1, Nx-1
    rhs_v(i) = c_wave**2 * (my_u(i+1) - 2.0d0*my_u(i) + my_u(i-1)) / dx**2
  end do

  ! Condiciones de frontera periódicas
  rhs_v(0)  = c_wave**2 * (my_u(1)    - 2.0d0*my_u(0)  + my_u(Nx-1)) / dx**2
  rhs_v(Nx) = c_wave**2 * (my_u(1)    - 2.0d0*my_u(Nx) + my_u(Nx-1)) / dx**2

end subroutine
