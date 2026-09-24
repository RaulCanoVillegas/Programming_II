program newtonraphson_ln
implicit none

real(8) :: x0, x1, tol, error
integer :: k

tol = 1.0d-8
x0 = 3.0d0
k = 0

open(10,file="newtonraphson_ln.dat")

do
    	x1 = x0 - (log(x0)-0.7d0)/(1.d0/x0)

    	if (k == 0) then
        	error = 1.0d0
    	else
        	error = abs((x1 - x0)/x0)
   	end if

    	write(10,*) k, x1, error

    	if(error < tol) exit

    	x0 = x1
    	k = k + 1
end do

print*, "Iteraciones:", k
print*, "Raiz:", x1

close(10)
end program
