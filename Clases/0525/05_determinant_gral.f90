! Code calculate determinant

program matrix_vector
implicit none

integer, parameter :: n = 5	! matrix size
integer	:: i
real(kind=8) :: A(n,n), det

! Define matrix 
A = reshape((/ &
	 1.0d0,  2.0d0,  0.0d0, -1.0d0,  3.0d0, &
         0.0d0,  1.0d0,  4.0d0,  2.0d0, -2.0d0, &
        -3.0d0,  0.0d0,  1.0d0,  5.0d0,  1.0d0, &
	 2.0d0, -1.0d0,  3.0d0,  0.0d0,  4.0d0, &
	 1.0d0,  2.0d0, -2.0d0,  1.0d0,  0.0d0 &
/), (/n,n/), order=(/2,1/))

det = determinant(A,n)

print *, "Determinant =", det

contains

! recursive
recursive function determinant(M,n) result(det)

implicit none
integer, intent(in) :: n
real(kind=8), intent(in) :: M(n,n)
real(kind=8) :: det, minor(n-1,n-1)
integer :: i,j,col,mi,mj

! Base cases
if (n == 1) then
	det = M(1,1)
elseif (n == 2) then
	det = M(1,1)*M(2,2) - M(1,2)*M(2,1)
else
	det = 0.0d0
	! expansion by cofactors
	do col = 1, n
		mi = 0
		do i = 2, n
			mi = mi + 1
			mj = 0
			do j = 1, n
				if (j /= col) then
					mj = mj + 1
					minor(mi,mj) = M(i,j)
				end if
			end do
		end do
		det =  det + (-1.0d0)**(1+col) * M(1,col) * determinant(minor,n-1)
	end do
end if

end function determinant
end program

