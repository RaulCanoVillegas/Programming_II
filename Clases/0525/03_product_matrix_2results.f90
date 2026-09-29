! Código de una matriz por un vector

program matrix_vector
implicit none

integer, parameter :: n = 5, m=5, p=5	! matrix size
integer	:: i, j, k
real(kind=8) :: A(n,m), B(m,p), AB(n,p), BA(m,m)

! Define matrix A 
A = reshape((/ &
	 1.0d0,  2.0d0,  0.0d0, -1.0d0,  3.0d0, &
         0.0d0,  1.0d0,  4.0d0,  2.0d0, -2.0d0, &
        -3.0d0,  0.0d0,  1.0d0,  5.0d0,  1.0d0, &
	 2.0d0, -1.0d0,  3.0d0,  0.0d0,  4.0d0, &
	 1.0d0,  2.0d0, -2.0d0,  1.0d0,  0.0d0 &
/), (/n,m/), order=(/2,1/))

! Define matrix B
B = reshape((/ &
	 1.0d0,  0.0d0,  2.0d0, -1.0d0,  3.0d0, &
         2.0d0,  1.0d0,  0.0d0,  4.0d0, -2.0d0, &
        -1.0d0,  3.0d0,  1.0d0,  0.0d0,  2.0d0, &
	 0.0d0, -2.0d0,  5.0d0,  1.0d0,  1.0d0, &
	 4.0d0,  1.0d0, -1.0d0,  2.0d0,  0.0d0 &
/), (/n,m/), order=(/2,1/))

! Matrix-matrix multiplication
AB = 0.0d0
do i = 1, n
	do j = 1, p
	        do k = 1, m
	           AB(i,j) = AB(i,j) + A(i,k)*B(k,j)
	        end do
	end do
end do

BA = 0.0d0
do i = 1, m
    do j = 1, m
        do k = 1, n
            BA(i,j) = BA(i,j) + B(i,k)*A(k,j)
        end do
    end do
end do

! Print
print *, "Matrix AB:"
do i = 1, n
        print *, AB(i,:)
end do

print *, "Matrix BA:"
do i = 1, m
        print *, BA(i,:)
end do

end program
