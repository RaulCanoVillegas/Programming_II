program numerical_derivate
! This program builds a 1D mesh

implicit none

real(kind=8), allocatable, dimension(:) :: x,f,df,error
integer :: i,N
real(kind=8) :: xmin, xmax, dx

! Initial conditions
xmin = -5.0
xmax =  5.0
N    =  800

! Allocate memory
allocate(x(0:N),f(0:N),df(0:N),error(0:N))

! Separation between points
dx = (xmax - xmin) / dble(N)

! Build x
do i=0,N
	x(i) = xmin + dble(i) * dx
end do

! Evaluate function
f = sin(x)

! Calculate derivate
do i=1,N-1
	df(i) = 0.5d0 * ( f(i+1) - f(i-1) ) / dx
end do
df(0) =-0.5d0 * ( f(2) - 4.0d0*f(1 ) + 3.0d0*f(0) ) / dx	! left border
df(N) = 0.5d0 * ( f(N-2) - 4.0d0*f(N-1) + 3.0d0*f(N) ) / dx	! right border
error = df - cos(x)

! Save
open(1,file='N800.dat')
do i=0,N
	write(1,*) x(i),f(i),df(i),error(i)
end do
close(1)

end program
