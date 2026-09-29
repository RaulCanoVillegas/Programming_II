! Code gauss_method 

program gauss_method
implicit none

integer, parameter :: n = 5
integer :: i, j, k
real(kind=8) :: A(n,n), b(n), Aug(n,n+1), x(n), factor, suma

! Define matrix
A = reshape((/ &
	 1.0d0,  2.0d0,  0.0d0, -1.0d0,  3.0d0, &
         0.0d0,  1.0d0,  4.0d0,  2.0d0, -2.0d0, &
        -3.0d0,  0.0d0,  1.0d0,  5.0d0,  1.0d0, &
	 2.0d0, -1.0d0,  3.0d0,  0.0d0,  4.0d0, &
	 1.0d0,  2.0d0, -2.0d0,  1.0d0,  0.0d0 &
/), (/n,n/), order=(/2,1/))

! Right-hand size vector
b = (/ 15.0d0, -6.0d0, -5.0d0, 12.0d0, 4.0d0 /) 
! (/ 13.0d0, 0.0d0, -1.0d0, 16.0d0, 4.0d0 /)
! (/ 4.0d0, 16.0d0, -1.0d0, 0.0d0, 13.0d0 /)


! Build argument matrix [A|b]
Aug(:,1:n) = A
Aug(:,n+1) = b

! Gaussian elimination
do k = 1, n-1
	if (abs(Aug(k,k)) < 1.0d-12) then
		print *, "Zero pivot encountered"
		stop
	end if
	do i = k+1, n
		factor = Aug(i,k)/Aug(k,k)
		do j = k, n+1
			Aug(i,j) = Aug(i,j) - factor * Aug(k,j)
		end do
	end do
end do

! Triangular matrix
x(n) = Aug(n,n+1)/Aug(n,n)
do i = n-1, 1, -1
	suma = 0.0d0
	do j = i+1, n
		suma = suma + Aug(i,j)*x(j)
	end do
	x(i) = (Aug(i,n+1) - suma)/Aug(i,i)
end do

! Print upper triangular matrix
print *, ""
print *, " Upper triangular matrix:"
do i = 1, n
	write(*,'(6F12.6)') Aug(i,:)
end do

! Back substitution
x(n) = Aug(n,n+1)/Aug(n,n)
do i = n-1, 1, -1
	suma = 0.0d0
	do j = i+1, n
		suma = suma + Aug(i,j)*x(j)
	end do
	x(i) = (Aug(i,n+1) - suma)/Aug(i,i)	!
end do

! Print solution
print *, ""
print *, "Solution:"
print *, ""
do i = 1, n
	print *, "x(",i,") = ", x(i)
end do

end program
