program biseccion
implicit none

real(8) :: raiz1, raiz2

print*, "Function: f(x) = x**2 - 2"
print*, "LOCATOR"
print*, " "
print*, "x1 [-2,0]"
call biseccion_metodo(-2.0d0, 0.0d0, raiz1)

print*, "x2 [0,2]"
call biseccion_metodo(0.0d0, 2.0d0, raiz2)

print*, " "
print*, "x1 =", raiz1
print*, "x2 =", raiz2

contains

!____________________________________________________
subroutine biseccion_metodo(a, b, raiz)
implicit none

real(8), intent(in)  :: a, b
real(8), intent(out) :: raiz
real(8) :: ai, bi, ci, fa, fb, fc
real(8) :: tol
integer :: iter, max_iter

tol = 1.0d-8
max_iter = 50

ai = a
bi = b

fa = f(ai)
fb = f(bi)

if (fa*fb > 0.0d0) then
    print*, "There isn't change of sign"
    raiz = 0.0d0
    return
end if

print*, "         iter         a                          b                         c                        f(c)"
print*, "--------------------------------------------------------------------------------------------------------------------------"

do iter = 1, max_iter

    ci = 0.5d0*(ai + bi)
    fc = f(ci)

    print*, iter, ai, bi, ci, fc

    if (abs(fc) < tol) exit

    if (fa*fc < 0.0d0) then
        bi = ci
        fb = fc
    else
        ai = ci
        fa = fc
    end if

end do

raiz = ci

end subroutine bisectiom_mothod
!____________________________________________________

! Function
function f(x)
implicit none
real(kind=8) :: f
real(kind=8), intent(in) :: x

f = x**2 - 2.0d0

end function f

end program 
