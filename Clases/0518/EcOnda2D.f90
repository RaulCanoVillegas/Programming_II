module numbers

integer :: Nx,Ny,Nt
real(kind=8) :: t,dt,tf,dx,dy,xmin,xmax,ymin,ymax
real(kind=8), allocatable, dimension(:,:) :: x,y
real(kind=8), allocatable, dimension(:,:) :: rhs_phi, phi, phi_p
real(kind=8), allocatable, dimension(:,:) :: rhs_psix,psix,psix_p
real(kind=8), allocatable, dimension(:,:) :: rhs_psiy,psiy,psiy_p
real(kind=8), allocatable, dimension(:,:) :: rhs_pi,  pi,  pi_p

end module

! -------------------

program wave

use numbers
implicit none

integer :: i,j,n

Nx   = 100
Ny   = 100
tf   = 2.0d0
xmin = -1.0d0
xmax =  1.0d0
ymin = -1.0d0
ymax =  1.0d0

allocate(x(0:Nx,0:Ny),y(0:Nx,0:Ny))
allocate(rhs_phi(0:Nx,0:Ny) ,phi(0:Nx,0:Ny) ,phi_p(0:Nx,0:Ny) )
allocate(rhs_psix(0:Nx,0:Ny),psix(0:Nx,0:Ny),psix_p(0:Nx,0:Ny))
allocate(rhs_psiy(0:Nx,0:Ny),psiy(0:Nx,0:Ny),psiy_p(0:Nx,0:Ny))
allocate(rhs_pi(0:Nx,0:Ny)  ,pi(0:Nx,0:Ny)  ,pi_p(0:Nx,0:Ny)  )

dx = ( xmax - xmin ) / dble(Nx)
dy = ( ymax - ymin ) / dble(Ny)
dt = 0.125 * dx
do i=0,Nx
  x(i,:) = xmin + dble(i) * dx
end do 
do j=0,Ny
  y(:,j) = ymin + dble(j) * dy
end do 
Nt = int(tf/dt)

open(1,file='EcOnda2D.dat')
  
! Initial conditions
t    = 0.0d0
phi  = exp( -(x**2+y**2) / 0.01 )
psix = -2.0*x*phi/0.01
psiy = -2.0*y*phi/0.01
pi   = 0.0

  do i=0,Nx
    do j=0,Ny
      write(1,*) x(i,j), y(i,j),phi(i,j)
    end do
    write(1,*)
  end do

  do n=1,Nt
    t      = t + dt
    phi_p  = phi
    psix_p = psix
    psiy_p = psiy
    pi_p   = pi

    call calcrhs( phi_p,psix_p,psiy_p,pi_p )	! Calculate dU/dt = F(U) - con U = (phi,psix,psiy,pi)
    phi  = phi_p  + rhs_phi  * dt    		! Predictor - Euler
    psix = psix_p + rhs_psix * dt    
    psiy = psiy_p + rhs_psiy * dt    
    pi   = pi_p   + rhs_pi   * dt

    call calcrhs( phi,psix,psiy,pi   ) 			! Evaluate RHS using the predictor
    phi  = 0.5d0 * ( phi_p  + phi  + rhs_phi  * dt )	! Corrector RK2 - U**(n+1) = 1/2 * (U**n + U**(*) + dtF(U**(*))
    psix = 0.5d0 * ( psix_p + psix + rhs_psix * dt )
    psiy = 0.5d0 * ( psiy_p + psiy + rhs_psiy * dt )
    pi   = 0.5d0 * ( pi_p   + pi   + rhs_pi   * dt )

    if ( mod(n,10).eq.0 ) then
      write(1,*)
      write(1,*)
      do i=0,Nx
        do j=0,Ny
          write(1,*) x(i,j),y(i,j), phi(i,j)
        end do
        write(1,*)
      end do
      print *, n,t
    end if
  end do
close(1)

end program

! ----------------------
subroutine calcrhs(my_phi,my_psix,my_psiy,my_pi)	! Calculate special derivates

use numbers
implicit none

real(kind=8), dimension(0:Nx,0:Ny), intent(in) :: my_phi,my_psix,my_psiy,my_pi
integer :: i,j

rhs_phi = my_pi

do i=1,Nx-1
  do j=1,Ny-1
    rhs_psix(i,j) = 0.5d0 * ( my_pi(i+1,j)   - my_pi(i-1,j)   ) / dx	! diferencia centrada de segundo origen
    rhs_psiy(i,j) = 0.5d0 * ( my_pi(i,j+1)   - my_pi(i,j-1)   ) / dy
    rhs_pi(i,j)   = 0.5d0 * ( my_psix(i+1,j) - my_psix(i-1,j) ) / dx &
                  + 0.5d0 * ( my_psiy(i,j+1) - my_psiy(i,j-1) ) / dy
  end do
end do

rhs_psix(0,:)  = 0.5d0 * ( my_pi(1,:) - my_pi(Nx-1,:) ) / dx	! Condiciones de fronteraz
rhs_psix(Nx,:) = 0.5d0 * ( my_pi(1,:) - my_pi(Nx-1,:) ) / dx
rhs_psiy(:,0)  = 0.5d0 * ( my_pi(:,1) - my_pi(:,Ny-1) ) / dy
rhs_psiy(:,Ny) = 0.5d0 * ( my_pi(:,1) - my_pi(:,Ny-1) ) / dy

do j=1,Ny-1
  rhs_pi(0,j)  = 0.5d0 * ( my_psix(1,j)    - my_psix(Nx-1,j) ) / dx &
               + 0.5d0 * ( my_psiy(0,j+1)  - my_psiy(0,j-1) )  / dy
  rhs_pi(Nx,j) = 0.5d0 * ( my_psix(1,j)    - my_psix(Nx-1,j) ) / dx &
               + 0.5d0 * ( my_psiy(Nx,j+1) - my_psiy(Nx,j-1) ) / dy
end do
do i=1,Nx-1
  rhs_pi(i,0)  = 0.5d0 * ( my_psix(i+1,0)  - my_psix(i-1,0)  ) / dx &
               + 0.5d0 * ( my_psiy(i,1)    - my_psiy(i,Ny-1) ) / dy
  rhs_pi(i,Ny) = 0.5d0 * ( my_psix(i+1,Ny) - my_psix(i-1,Ny) ) / dx &
               + 0.5d0 * ( my_psiy(i,1)    - my_psiy(i,Ny-1) ) / dy
end do
rhs_pi(0,0)   = 0.5d0 * ( my_psix(1 ,0 ) - my_psix(Nx-1,0   ) ) / dx &
              + 0.5d0 * ( my_psiy(0 ,1 ) - my_psiy(0   ,Ny-1) ) / dy
rhs_pi(0,Ny)  = 0.5d0 * ( my_psix(1 ,Ny) - my_psix(Nx-1,Ny  ) ) / dx &
              + 0.5d0 * ( my_psiy(0 ,1 ) - my_psiy(0   ,Ny-1) ) / dy
rhs_pi(Nx,0)  = 0.5d0 * ( my_psix(1 ,0 ) - my_psix(Nx-1,0   ) ) / dx &
              + 0.5d0 * ( my_psiy(Nx,1 ) - my_psiy(Nx  ,Ny-1) ) / dy
rhs_pi(Nx,Ny) = 0.5d0 * ( my_psix(1 ,Ny) - my_psix(Nx-1,Ny  ) ) / dx &
              + 0.5d0 * ( my_psiy(Nx,1 ) - my_psiy(Nx  ,Ny-1) ) / dy



end subroutine





























