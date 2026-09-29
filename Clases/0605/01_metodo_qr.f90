	! Code qr_method

program qr_method

implicit none
integer, parameter :: n = 4
real :: A(n,n), Q(n,n), R(n,n), Temp(n,n), tol, error
integer :: i, j, k, iter, max_iter

A(1,:) = (/ 4.0d0,  1.0d0,  0.0d0, 0.0d0 /)
A(2,:) = (/ 1.0d0,  4.0d0,  1.0d0, 0.0d0 /)
A(3,:) = (/ 0.0d0,  1.0d0,  4.0d0, 1.0d0 /)
A(4,:) = (/ 0.0d0,  0.0d0,  1.0d0, 4.0d0 /)

print *, "--- Matriz Original A ---"
do i = 1, n
        print '(10F8.3)', A(i,:)
end do
print *

! Parameters
max_iter = 200
tol = 1.0e-5

! QR METHOD
 do iter = 1, max_iter
! --- QR descomposition (Gram-Schmidt Modificade) ---
	Q = 0.0
        R = 0.0
        do j = 1, N
            	Q(:,j) = A(:,j)
            	do i = 1, j-1
                	R(i,j) = dot_product(Q(:,i), A(:,j))
                	Q(:,j) = Q(:,j) - R(i,j) * Q(:,i)
            	end do
            	R(j,j) = norm2(Q(:,j))
            	if (R(j,j) > 1.0e-7) then
                	Q(:,j) = Q(:,j) / R(j,j)
            	end if
        end do

        A = matmul(R, Q)

        error = 0.0
        do j = 1, n
            do i = j+1, n
                error = error + abs(A(i,j))
            end do
        end do

        if (error < tol) exit
end do

! Results
print '(A, I4, A)', "El método QR convergió en ", iter, " iteraciones."
print *
print *, "--- Matriz Diagonalizada (Autovalores en la diagonal) ---"
do i = 1, n
	print '(10F8.4)', A(i,:)
end do
print *

print *, "Los autovalores calculados son:"
do i = 1, n
	print '(A, I2, A, F8.4)', " l_", i, " = ", A(i,i)
end do

contains

real function norm2(vec)
real, intent(in) :: vec(:)
norm2 = sqrt(sum(vec**2))
end function norm2

end program 
