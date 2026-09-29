! Code gaussjordan_method 

program gauss_method
implicit none

integer, parameter :: n = 5
integer :: i, j, k
real(kind=8) :: A(n,n), b(n), Aug(n,n+1), x(n), factor, suma, pivot

! Define matrix
A = reshape((/ &
	 1.0d0,  2.0d0,  0.0d0, -1.0d0,  3.0d0, &
         0.0d0,  1.0d0,  4.0d0,  2.0d0, -2.0d0, &
        -3.0d0,  0.0d0,  1.0d0,  5.0d0,  1.0d0, &
	 2.0d0, -1.0d0,  3.0d0,  0.0d0,  4.0d0, &
	 1.0d0,  2.0d0, -2.0d0,  1.0d0,  0.0d0 &
/), (/n,n/), order=(/2,1/))

! Right-hand size vector
b = (/ 13.0d0, 0.0d0, -1.0d0, 16.0d0, 4.0d0 /) 
! (/ 15.0d0, -6.0d0, -5.0d0, 12.0d0, 4.0d0 /) 
! (/ 13.0d0, 0.0d0, -1.0d0, 16.0d0, 4.0d0 /)


! Build argument matrix [A|b]
Aug(:,1:n) = A
Aug(:,n+1) = b

! Gaussia-Jordan elimination
do k = 1, n
	pivot = Aug(k,k)
	if (abs(pivot) < 1.0d-12) then
		print *, "Zero pivot encountered"
		stop
	end if
	! normalize pivot row
	do j = 1, n+1
		Aug(k,j) = Aug(k,j)/pivot
	end do
	! eliminate above and below
	do i = 1, n
		if (i /= k) then
			factor = Aug(i,k)
			do j = 1, n+1
				Aug(i,j) = Aug(i,j) - factor*Aug(k,j)
			end do
		end if
	end do
end do

! Print reduced matrix
print *, ""
print *, " Reduced matrix:"
do i = 1, n
	write(*,'(6F12.6)') Aug(i,:)
end do

! Extract solution
x = Aug(:,n+1)

! Print solution
print *, ""
print *, "Solution:"
print *, ""
do i = 1, n
	print *, "x(",i,") = ", x(i)
end do

end program
