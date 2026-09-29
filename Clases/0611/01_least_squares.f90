! Code least_squares_qr

program least_squares_qr

implicit none
integer, parameter :: m = 6   ! number of equations (data points)
integer, parameter :: n = 2   ! number of unknowns (e.g. 2 for a line y = a + b*x)
integer :: i, j
real(kind=8) :: A(m,n), b(m), Q(m,n), R(n,n), Qtb(n), x(n)

! ---- Example data: fit y = a + b*x ----
! x values: 1,2,3,4,5,6 ; y values: 2.1, 2.9, 3.6, 4.2, 5.1, 5.8
A(:,1) = (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)   ! column of 1's (intercept)
A(:,2) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)   ! x values

b = (/ 2.1d0, 2.9d0, 3.6d0, 4.2d0, 5.1d0, 5.8d0 /)

print *, "Matrix A:"
do i = 1, m
	print '(10F8.3)', A(i,:)
end do
print *

print *, "Vector b:"
do i = 1, m
	print '(F8.3)', b(i)
end do
print *

! ---- QR factorization (modified Gram-Schmidt) ----
Q = 0.0d0
R = 0.0d0
do j = 1, n
	Q(:,j) = A(:,j)
	do i = 1, j-1
		R(i,j) = dot_product(Q(:,i), A(:,j))
		Q(:,j) = Q(:,j) - R(i,j) * Q(:,i)
	end do
	R(j,j) = norm2_f(Q(:,j))
	if (R(j,j) > 1.0d-7) Q(:,j) = Q(:,j) / R(j,j)
end do

print *, "Matrix Q:"
do i = 1, m
	print '(10F8.4)', Q(i,:)
end do
print *

print *, "Matrix R:"
do i = 1, n
	print '(10F8.4)', R(i,:)
end do
print *

! ---- Compute Q^T * b ----
Qtb = matmul(transpose(Q), b)

print *, "Q^T b:"
do i = 1, n
	print '(F8.4)', Qtb(i)
end do
print *

! ---- Solve R x = Q^T b by back substitution ----
do i = n, 1, -1
	x(i) = Qtb(i)
	do j = i+1, n
		x(i) = x(i) - R(i,j) * x(j)
	end do
	if (abs(R(i,i)) > 1.0d-10) then
		x(i) = x(i) / R(i,i)
	else
		print *, "WARNING: R has a (near) zero diagonal element, system may be singular"
	end if
end do

print *, "Least squares solution x:"
do i = 1, n
	print '(A, I2, A, F8.4)', "  x_", i, " = ", x(i)
end do
print *

! ---- Residuals ----
print *, "Fitted values (A*x) and residuals (b - A*x):"
do i = 1, m
	print '(A, F8.4, A, F8.4)', "  fitted = ", dot_product(A(i,:), x), &
	       "   residual = ", b(i) - dot_product(A(i,:), x)
end do

contains

! ---- Euclidian norm ----
real(kind=8) function norm2_f(vec)
	real(kind=8), intent(in) :: vec(:)
	norm2_f = sqrt(sum(vec**2))
end function norm2_f

end program
