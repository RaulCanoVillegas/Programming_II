program localizador_raices
implicit none

real :: a,b,dx,x
real :: f1,f2,raiz
real :: tol
real :: a_int,b_int
integer :: iter,maxit
logical :: raiz_doble_encontrada

a = -2.0
b = 6.0
dx = 0.2
tol = 1e-10
maxit = 100
x = a
raiz_doble_encontrada = .false.

print *, "Localizacion de raices de f(x)=(x-2)^4-10(x-2)^2"
print *, " "

do while (x < b)
    	f1 = f(x)
    	f2 = f(x+dx)
    	! raices simples
    	if (f1*f2 < 0.0) then
        	a_int = x
        	b_int = x + dx
		call biseccion(a_int,b_int,tol,maxit,raiz,iter)
		print *, "Raiz encontrada:", raiz
        	print *, "Iteraciones:", iter
        	print *, " "
    	end if
	! raiz doble cerca de x=2
	if (abs(f1) < 1e-3 .and. .not. raiz_doble_encontrada) then
		call newton_modificado(x,tol,maxit,raiz,iter)
		print *, "Raiz doble encontrada en:", raiz
        	print *, "Iteraciones:", iter
        	print *, " "
        	raiz_doble_encontrada = .true.
    	end if
	x = x + dx
end do


contains

real function f(x)
real :: x
f = (x-2.0)**4 - 10.0*(x-2.0)**2
end function

real function df(x)
real :: x
df = 4.0*(x-2.0)**3 - 20.0*(x-2.0)
end function

real function ddf(x)
real :: x
ddf = 12.0*(x-2.0)**2 - 20.0
end function

subroutine biseccion(a,b,tol,maxit,raiz,iter)
real :: a,b,c,tol,raiz
integer :: iter,maxit
iter = 0
do while (abs(b-a) > tol .and. iter < maxit)
    	c = (a+b)/2.0
    	if (f(a)*f(c) < 0.0) then
        	b = c
    	else
        	a = c
    	end if
    	iter = iter + 1
end do
raiz = (a+b)/2.0
end subroutine

subroutine newton_modificado(x0,tol,maxit,raiz,iter)
real :: x0,tol,raiz
integer :: iter,maxit
real :: x,x_new
real :: fx,dfx,ddfx,denom,dxn
x = x0
do iter = 1,maxit
    	fx = f(x)
    	dfx = df(x)
    	ddfx = ddf(x)
    	denom = dfx**2 - fx*ddfx
	if (abs(denom) < 1e-12) exit
    	dxn = (fx*dfx)/denom
	x_new = x - dxn
	if (abs(x_new-x) < tol) exit
	x = x_new
end do
raiz = x
end subroutine

end program localizador_raices
