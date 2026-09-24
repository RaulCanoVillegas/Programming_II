program simple_harmonic_oscillator

implicit none
integer :: i,N
real(kind=8) :: t,dt,k,m,x0,y0,x_p,x,y,y_p,xE,error,tf
real(kind=8) :: k1x,k2x,k1y,k2y

N  = 100
tf = 25.0
k  = 2.0
m  = 0.5
x0 = 1.0
y0 = 0.0
dt = tf / dble(N)

t  = 0.0
x  = x0
y  = y0
xE = x0

error = x-xE
print *, t,x,y,xE,error

do i=0,N
  t   = t + dt
  x_p = x
  y_p = y
  ! --------------------------------------------------------
  k1x = y_p
  k1y = -k/m * x_p
  ! temporal euler
  x   = k1x * dt + x_p	
  y   = k1y * dt + y_p
  ! new slopes
  k2x = y
  k2y = -k/m * x
  x   = k2x * dt + x_p
  y   = k2y * dt + y_p
  ! --------------------------------------------------------
  xE = cos( 2.0 * t )
  error = x - xE
  print*, t,x,y,xE,error
end do

end program
