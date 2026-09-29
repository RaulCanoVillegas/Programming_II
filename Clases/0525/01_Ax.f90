! Código de una matriz por un vector

program matrix_vector
implicit none

integer, parameter :: n = 5, m=5	! matrix size
integer	:: i, j
real(kind=8) :: A(n,m), x(n), b(n)

! Define matrix 
A = reshape((/ &
	 1.0d0,  2.0d0,  0.0d0, -1.0d0,  3.0d0, &
         0.0d0,  1.0d0,  4.0d0,  2.0d0, -2.0d0, &
        -3.0d0,  0.0d0,  1.0d0,  5.0d0,  1.0d0, &
	 2.0d0, -1.0d0,  3.0d0,  0.0d0,  4.0d0, &
	 1.0d0,  2.0d0, -2.0d0,  1.0d0,  0.0d0 &
/), (/n,m/), order=(/2,1/))

! Define vector
x = (/ 1.0d0, -1.0d0, 2.0d0, 0.0d0, 3.0d0 /)

b = 0.0d0

! Matrix-vector multiplication
do i = 1, n
	do j = 1, m
            	b(i) = b(i) + A(i,j) * x(j)
        end do
end do

! Print
print *, "Resultado:"
    	do i = 1, n
        	print *, b(i)
    	end do

end program
