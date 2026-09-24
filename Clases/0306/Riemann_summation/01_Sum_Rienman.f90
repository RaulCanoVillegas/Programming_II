program Sum_Rienman
implicit none

integer :: i,N
real(kind=8) :: xmin,xmax,dx,Integral
real(kind=8), allocatable :: x(:),f(:),df(:),error(:)

xmin = -5.0d0
xmax = 5.0d0
N    = 100
allocate(x(0:N),f(0:N),df(0:N),error(0:N))
dx = (xmax - xmin)/dble(N)	! discretization

do i = 0,N
    x(i) = xmin + dble(i)*dx	! Build points
end do

! Function
f = exp(-x**2)

do i = 1,N-1
    df(i) = 0.5d0 * ( f(i+1) - f(i-1) ) / dx	! numerical_derivate
end do

df(0) = (-3.0d0*f(0) + 4.0d0*f(1) - f(2)) / (2.0d0*dx)		! border_x0
df(N) = (3.0d0*f(N) - 4.0d0*f(N-1) + f(N-2)) / (2.0d0*dx)	! border_xN
error = df + 2.0d0 * x * f

Integral = 0.0d0
do i = 0,N-1
    	Integral = Integral + f(i)*dx
end do

print*, 'Approximate integral=', Integral
print*, 'Exact integral =', sqrt(acos(-1.0d0))
print*, 'Error =', Integral - sqrt(acos(-1.0d0))

open(1,file='Sum_Riemann.dat',status='replace')
do i = 0,N
    write(1,'(5F15.8)') x(i), f(i), df(i), -2.0d0*x(i)*f(i), error(i)
end do
close(1)

end program 
