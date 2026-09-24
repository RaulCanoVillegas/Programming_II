program raices_multiples
implicit none

real :: x0, raiz
real :: tol
integer :: iter, maxit
integer :: m_est

x0 = 1.5
tol = 1e-10
maxit = 100

call newton_modificado(x0,tol,maxit,raiz,iter)
m_est = nint(estimar_multiplicidad(raiz))

print *, "Raiz multiple encontrada:", raiz
print *, "Iteraciones:", iter
print *, "Multiplicidad aproximada:", m_est


contains

real function f(x)
real :: x
f = (x-2.0)**3
end function

real function df(x)
real :: x
df = 3.0*(x-2.0)**2
end function

real function ddf(x)
real :: x
ddf = 6.0*(x-2.0)
end function


subroutine newton_modificado(x0,tol,maxit,raiz,iter)

real :: x0,tol,raiz
integer :: iter,maxit
real :: x,x_new
real :: fx,dfx,ddfx,denom,dx

x = x0

do iter = 1,maxit
    	fx = f(x)
    	dfx = df(x)
    	ddfx = ddf(x)
	denom = dfx**2 - fx*ddfx
	if (abs(denom) < 1e-12) exit
	dx = (fx*dfx)/denom
    	x_new = x - dx
    	if (abs(x_new-x) < tol) exit
    	x = x_new
end do
raiz = x
end subroutine

real function estimar_multiplicidad(x)
real :: x
real :: fx,dfx,ddfx,R

fx = f(x)
dfx = df(x)
ddfx = ddf(x)

if (abs(dfx) < 1e-12) then
    	estimar_multiplicidad = 3.0
else
    	R = (fx*ddfx)/(dfx**2)
    	estimar_multiplicidad = 1.0/(1.0 - R)
end if
end function

end program 
