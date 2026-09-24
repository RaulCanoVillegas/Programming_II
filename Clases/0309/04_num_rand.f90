program num_aleatorio
implicit none

real(kind=8)                            :: u, dx
real(kind=8), allocatable, dimension(:) :: x 
integer, allocatable, dimension(:)      :: histograma
integer                                 :: i,j,N,N_cells

N       = 1000000
N_cells = 1000
dx      = 1.0d0 / dble(N_cells)

allocate(x(0:N_cells),histograma(0:N_cells-1))

do i=0,N_cells
	x(i) = dble(i) * dx
end do

histograma = 0
do j=0,N
	call random_number(u)
	do i=0,N_cells-1
		if ( ( u.ge.x(1) ).and.( u.lt.x(i+1) ) ) then
			histograma(i) = histograma(i) + 1
			exit
		end if
	end do
end do

open(1,file='Num_rand.dat')
do i=0,N_cells-1
	write(1,*) x(i)+dx/2, histograma(i)
end do
close(1)

end program
