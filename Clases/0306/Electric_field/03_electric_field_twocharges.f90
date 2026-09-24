program electric_field_anim
  implicit none

  integer :: i, j, k, Nx, Ny, Nt
  real(kind=8) :: xmin, xmax, ymin, ymax, dx, dy
  real(kind=8) :: t, dt, pi
  real(kind=8) :: r, escala
  real(kind=8), allocatable :: x(:,:), y(:,:)
  real(kind=8), allocatable :: vx(:,:), vy(:,:)
  real(kind=8) :: xq, yq   ! position of the charge
  real(kind=8) :: xq2, yq2 ! position of the second charge (opposite)
  real(kind=8) :: ex, ey   ! total field

  ! Parametes
  Nx     = 25
  Ny     = 25
  Nt     = 72              
  xmin   = -2.0d0
  xmax   =  2.0d0
  ymin   = -2.0d0
  ymax   =  2.0d0
  pi     = acos(-1.0d0)
  escala = 0.3d0           
  
  ! Separation between points
  dx = (xmax - xmin) / dble(Nx)
  dy = (ymax - ymin) / dble(Ny)
  dt = 2.0d0 * pi / dble(Nt)

  ! dynamics arrays
  allocate(x(0:Nx, 0:Ny))
  allocate(y(0:Nx, 0:Ny))
  allocate(vx(0:Nx, 0:Ny))
  allocate(vy(0:Nx, 0:Ny))

  ! Built mesh
  do i = 0, Nx
    do j = 0, Ny
      x(i,j) = xmin + i * dx
      y(i,j) = ymin + j * dy
    end do
  end do

  ! save
  open(10, file='animated_field.dat', status='replace')

  do k = 0, Nt - 1
    t = k * dt

    ! Charge +1 orbits counterclockwise to radio 0.6
    xq  =  0.6d0 * cos(t)
    yq  =  0.6d0 * sin(t)

    ! Charge -1 diametrically opposite (rotating dipole)
    xq2 = -0.6d0 * cos(t)
    yq2 = -0.6d0 * sin(t)

    do i = 0, Nx
      do j = 0, Ny

        ! Charge's field +1
        r = sqrt((x(i,j) - xq)**2 + (y(i,j) - yq)**2)
        if (r > 0.12d0) then
          ex =  escala * (x(i,j) - xq) / r**3
          ey =  escala * (y(i,j) - yq) / r**3
        else
          ex = 0.0d0
          ey = 0.0d0
        end if

        ! Sum charge's field -1
        r = sqrt((x(i,j) - xq2)**2 + (y(i,j) - yq2)**2)
        if (r > 0.12d0) then
          ex = ex - escala * (x(i,j) - xq2) / r**3
          ey = ey - escala * (y(i,j) - yq2) / r**3
        end if

        ! Saturate magnitude so the arrows aren't gigantic
        r = sqrt(ex**2 + ey**2)
        if (r > 0.6d0) then
          ex = ex / r * 0.6d0
          ey = ey / r * 0.6d0
        end if

        vx(i,j) = ex
        vy(i,j) = ey

        write(10,*) x(i,j), y(i,j), vx(i,j), vy(i,j)
      end do
      write(10,*)   
    end do
    
    write(10,*)     
    write(10,*)

  end do

  close(10)

  deallocate(x, y, vx, vy)
end program electric_field_anim
