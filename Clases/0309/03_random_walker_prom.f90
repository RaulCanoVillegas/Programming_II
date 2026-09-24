program random_walker_prom
implicit none

integer                                 :: i,j,steps,walkers
real(kind=8)                            :: f,x
real(kind=8), allocatable, dimension(:) :: avg

steps   = 10000
walkers = 1000
allocate(avg(0:steps))
avg = 0.0d0

open(1,file='av_random_walker.dat')
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
		avg(i) = avg(i) + f
	end do
	write(1,*)
	write(1,*)
end do
close(1)
avg = avg / dble(walkers+1)

open(2,file='prom')
do i=0,steps
	write(2,*) i, avg(i)
end do
close(2)

! plot "salida4" w l lc rgb "purple" notitle, \
    ! "prom" w l lw 3 lc title "Promedio"

end program
