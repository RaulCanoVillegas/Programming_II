program interpolacion_sin
implicit none

real(kind=8) :: xt, xu, xr, xr_old, tol, error
integer :: k

tol = 1.0d-8

xt = 0.5d0
xu = 1.0d0
xr_old = 0.0d0
k = 0

open(10,file="interpolacion_sin.dat")

do
    	xr = xu - ((sin(xu)-xu**2)*(xt-xu)) / &
         	((sin(xt)-xt**2)-(sin(xu)-xu**2))

    	if (k == 0) then
        	error = 1.0d0
    	else
        	error = abs((xr - xr_old)/xr_old)
    	end if

    	write(10,*) k, xr, error

    	if(error < tol) exit

    	if((sin(xt)-xt**2)*(sin(xr)-xr**2) < 0.d0) then
        	xu = xr
    	else
        	xt = xr
    	end if

    	xr_old = xr
    	k = k + 1
end do

print*, "Iteraciones:", k
print*, "Raiz:", xr

close(10)
end program

