program random_walker
implicit none

integer      :: i,j,steps,walkers
real(kind=8) :: f,x

steps   = 1000
walkers = 10

open(1,file='random_walker_simple.dat')
do j=0,walkers
	f = 0.0d0
	do i=0,steps
		call random_number(x)
		if ( x .lt. 0.5 ) then
			f = f + 1
		else if ( x .ge. 0.5) then
			f = f - 1
		end if
		write(1,*) i, f
	end do
	write(1,*)
	write(1,*)
end do
close(1)

end program
