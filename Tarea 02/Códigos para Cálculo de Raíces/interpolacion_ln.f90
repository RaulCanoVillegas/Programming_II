program interpolacion_ln
implicit none

real(kind=8) :: xt, xu, xr, xr_old, tol, error
integer :: k

tol = 1.0d-8
xt = 1.0d-6
xu = 3.0d0
xr_old = 0.0d0
k = 0

open(10,file="interpolacion_ln.dat")

do
    xr = xu - ((log(xu)-0.7d0)*(xt-xu)) / &
         ((log(xt)-0.7d0)-(log(xu)-0.7d0))

    if (k == 0) then
        error = 1.0d0
    else
        error = abs((xr - xr_old)/xr_old)
    end if

    write(10,*) k, xr, error

    if(error < tol) exit

    if((log(xt)-0.7d0)*(log(xr)-0.7d0) < 0.d0) then
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
