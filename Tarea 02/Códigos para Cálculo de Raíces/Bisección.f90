program biseccion
implicit none

real(kind=8) :: xt, xu, xr, xr_old, tol, error
integer :: k

tol = 1.0d-10
xt = 0.0d0
xu = 3.0d0
xr_old = 1.0d0
k = 0

open(10,file="biseccion.dat")

do
    	xr = (xt + xu)/2.0d0

    	if (k == 0) then
        	error = 1.0d0
    	else
        	error = abs((xr - xr_old)/xr_old)
    	end if

    	write(10,*) k, xr, error

    	if(error < tol) exit

    	if((xt**10 - 1.0d0)*(xr**10 - 1.0d0) < 0.0d0) then 
        	xu = xr
    	else
        	xt = xr
    	end if
    
    	xr_old = xr
    	k = k + 1
end do

print*, "Iteraciones:", k
print*, "Raiz aproximada:", xr

close(10)

end program
