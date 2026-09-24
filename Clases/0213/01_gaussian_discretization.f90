program Gaussian_discretization

implicit none

real(kind=8), allocatable, dimension(:) :: x,f
integer :: i,N
real(kind=8) :: xmin, xmax, dx

! Initial conditions
xmin = -5.0
xmax = 5.0
N = 100
allocate(x(0:N),f(0:N))
dx = ( xmax - xmin ) / dble(N)

! Build (variables)
do i=0,N
  x(i) = xmin + dble(i) * dx
end do

! Function
f = exp(-x**2)

! Save
open(1,file='gaussian_discretization.dat')
do i=0,N
  write(1,*) x(i),f(i)
end do
close(1)

end program
