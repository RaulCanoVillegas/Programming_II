! Code inverse_matrix

program inverse_matrix
implicit none

integer, parameter :: n = 5	! matrix size
integer	:: i, j, k, row, col, mi, mj
real(kind=8) :: A(n,n), C(n,n), Ainv(n,n), detA, minor(n-1,n-1), AAinv(n,n)

! Define matrix 
A = reshape((/ &
	 1.0d0,  2.0d0,  0.0d0, -1.0d0,  3.0d0, &
         0.0d0,  1.0d0,  4.0d0,  2.0d0, -2.0d0, &
        -3.0d0,  0.0d0,  1.0d0,  5.0d0,  1.0d0, &
	 2.0d0, -1.0d0,  3.0d0,  0.0d0,  4.0d0, &
	 1.0d0,  2.0d0, -2.0d0,  1.0d0,  0.0d0 &
/), (/n,n/), order=(/2,1/))

! Determinant
detA = determinant(A,n)
print *, "Determinant =", detA

! Cofactor_matrix
do row = 1, n 
	do col = 1, n
		mi = 0
		do i = 1, n
			if(i/=row) then
				mi = mi + 1
				mj = 0
				do j = 1, n
					if(j/=col) then
						mj = mj + 1
						minor(mi,mj) = A(i,j)
					end if
				end do
			end if
		end do
		C(row,col) = (-1.0d0)**(row+col) * determinant(minor,n-1)
	end do
end do

! Inverse_matrxi
Ainv = transpose(C) / detA

! Print
print *, ""
print *, "Inverse matrix"
print *, ""
do i = 1, n
	print *, Ainv(i,:)
end do

! Verification A*Ainv
AAinv = 0.0d0
do i = 1, n
 	do j = 1, n
		do k = 1, n
			AAinv(i,j) = AAinv(i,j) + Ainv(i,k)*A(k,j)
		end do
	end do
end do

print *, "A^{-1}A:"
do i = 1, n
	write(*,'(5ES15.5)') AAinv(i,:)
end do


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

