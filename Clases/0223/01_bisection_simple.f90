program bisection 
! It must be remembered that the bisection method requires a change of sign

implicit none
real(kind=8) :: raiz1, raiz2

! Interval
call biseccion_metodo(-2.0d0, 0.0d0, raiz1)
call biseccion_metodo(0.0d0, 2.0d0, raiz2)

print*, "x1 =", raiz1
print*, "x2 =", raiz2

contains

!______________________________________________________
subroutine biseccion_metodo(a, b, raiz)
implicit none

real(kind=8), intent(in)  :: a, b
real(kind=8), intent(out) :: raiz
real(kind=8) :: ai, bi, ci, fa, fb, fc
real(kind=8) :: tol
integer :: iter, max_iter

tol = 1.0d-8
max_iter = 100

! Copies
ai = a
bi = b

! Initial evalutation
fa = f(ai)
fb = f(bi)

! Verification (change sign)
if (fa*fb > 0.0d0) then
    print*, "There isn't sign change in the interval"
    raiz = 0.0d0
    return
end if

! Iterations
do iter = 1, max_iter

    ci = 0.5d0*(ai + bi)	! center point
    fc = f(ci)			! Evaluate function

    if (abs(fc) < tol) exit	! Stop condition

    if (fa*fc < 0.0d0) then	! Choose subinterval
        bi = ci
        fb = fc
    else
        ai = ci
        fa = fc
    end if

end do

raiz = ci

end subroutine biseccion_metodo
!______________________________________________________

function f(x)
implicit none
real(kind=8) :: f
real(kind=8), intent(in) :: x

! Function
f = x**2 - 2.0d0

end function f

end program 
