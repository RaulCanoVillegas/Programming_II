program montecarlo
implicit none

integer      :: i,N
real(kind=8) :: x, f, a, b, funtion, int_

a = 2
b = 3
N = 1000000

do i=1,N
	call random_number(x)
	x = ( b - a ) * x + a		! reescaldado
	funtion = x**2			
	int_ = int_ + (b-a) * funtion	! acumulation
end do

int_ = int_ / dble(N)			! average
print *, 'Integral = ', int_
print *, 'Exact integral = ', 9.0d0-8.0d0/3.0d0


end program
