! Code diagonalized_matrix

program diagonalized_matrix

implicit none
integer, parameter :: n = 4
integer :: i, j, iter, max_iter
real(kind=8) :: A(n,n), A_orig(n,n), Q(n,n), R(n,n), P(n,n), P_inv(n,n), D(n,n), eigenvalues(n), tol, error

A(1,:) = (/ 4.0d0, 1.0d0, 0.0d0, 0.0d0 /)
A(2,:) = (/ 1.0d0, 4.0d0, 1.0d0, 0.0d0 /)
A(3,:) = (/ 0.0d0, 1.0d0, 4.0d0, 1.0d0 /)
A(4,:) = (/ 0.0d0, 0.0d0, 1.0d0, 4.0d0 /)

A_orig = A 
print *, "Original matrix:"
do i = 1, n
	print '(10F8.3)', A(i,:)
end do
print *

! QR method
max_iter = 200
tol = 1.0d-5

do iter = 1, max_iter
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
        A = matmul(R, Q)

        error = 0.0d0
        do j = 1, n
            	do i = j+1, n
                	error = error + abs(A(i,j))
            	end do
        end do
        if (error < tol) exit
end do

print '(A, I4, A)', "Iteration", iter, "."
print *
print *, "Diagonalized matrix:"
do i = 1, n
	print '(10F8.4)', A(i,:)
end do
	print *

! eigenvalues
do i = 1, n
        eigenvalues(i) = A(i,i)
end do

print *, "Eigenvalues:"
    	do i = 1, n
        	print '(A, I2, A, F8.4)', "  l_", i, " = ", eigenvalues(i)
    	end do
    	print *

! Eigenvectors
call calc_eigenvectors(A_orig, eigenvalues, P, n)

print *, "P matrix (eigenvectors by columns):"
do i = 1, n
        print '(10F8.4)', P(i,:)
end do
print *

! P inverse
call gauss_jordan_inv(P, P_inv, n)

print *, "Matriz P_inv:"
do i = 1, n
	print '(10F8.4)', P_inv(i,:)
end do
print *

! D = P^{-1} A P 
D = matmul(P_inv, matmul(A_orig, P))

print *, "D = P^{-1} A P:"
do i = 1, n
        print '(10F8.4)', D(i,:)
end do

contains

! Euclidian norm
real(kind=8) function norm2_f(vec)
	real(kind=8), intent(in) :: vec(:)
        norm2_f = sqrt(sum(vec**2))
end function norm2_f

! Eigenvectors by hogeneous system
subroutine calc_eigenvectors(Amat, evals, Pmat, m)
        
integer, intent(in) :: m
real(kind=8), intent(in)    :: Amat(m,m), evals(m)
real(kind=8), intent(out)   :: Pmat(m,m)

real(kind=8) :: B(m,m), sys(m-1, m-1), rhs(m-1), factor, norma
integer :: k, i, j, p

do k = 1, m
! B = A - lambda_k * I
B = Amat
do i = 1, m
	B(i,i) = B(i,i) - evals(k)
end do

! Fix x_n = 1, build reduced system (m-1) x (m-1)
do i = 1, m-1
	do j = 1, m-1
		sys(i,j) = B(i,j)
        end do
        rhs(i) = -B(i,m)   ! Column m → RHS
end do

! Gaussian elimination with partial pivot
do p = 1, m-2
	! search pivot row
	do i = p+1, m-1
        	if (abs(sys(i,p)) > abs(sys(p,p))) then
                	sys( [p,i], :) = sys( [i,p], :)   ! swap rows
                        rhs( [p,i] ) = rhs( [i,p] )
                end if
	end do
        do i = p+1, m-1
		if (abs(sys(p,p)) < 1.0e-10) cycle
                    	factor = sys(i,p) / sys(p,p)
                    	sys(i,:) = sys(i,:) - factor * sys(p,:)
                    	rhs(i)   = rhs(i)   - factor * rhs(p)
                end do
	end do

	! Backward sustitution
        do i = m-1, 1, -1
		Pmat(i,k) = rhs(i)
                do j = i+1, m-1
                	Pmat(i,k) = Pmat(i,k) - sys(i,j) * Pmat(j,k)
                end do
                if (abs(sys(i,i)) > 1.0e-10) then
                    	Pmat(i,k) = Pmat(i,k) / sys(i,i)
                end if
	end do
        Pmat(m,k) = 1.0d0	! x_n	
end do

end subroutine calc_eigenvectors

! Inverse by Gauss-Jordan
subroutine gauss_jordan_inv(Mat, Inv, m)

integer, intent(in)  :: m
real(kind=8), intent(in)     :: Mat(m,m)
real(kind=8), intent(out)    :: Inv(m,m)

real(kind=8) :: Aug(m, 2*m), factor
integer :: i, j, p

! Matrix pĺus [Mat | I]
do i = 1, m
	Aug(i, 1:m)   = Mat(i,:)
        Aug(i, m+1:2*m) = 0.0d0
        Aug(i, m+i)   = 1.0d0
end do

! Forward elimination
do p = 1, m
! Partial pivot
	do i = p+1, m
		if (abs(Aug(i,p)) > abs(Aug(p,p))) then
			Aug([p,i],:) = Aug([i,p],:)
        	end if
	end do
	if (abs(Aug(p,p)) < 1.0d-10) then
		print *, "WARNING: Singular or near-singular matrix"
        	Inv = 0.0d0
        	return
	end if
	Aug(p,:) = Aug(p,:) / Aug(p,p)
	do i = 1, m
        	if (i == p) cycle
                	factor = Aug(i,p)
                	Aug(i,:) = Aug(i,:) - factor * Aug(p,:)
        end do
end do

Inv = Aug(1:m, m+1:2*m)
end subroutine gauss_jordan_inv

end program 
