program Newtons_law_of_cooling
implicit none

integer :: k,N
real(kind=8) :: t, dt, TempA, Temp0, Temp, TempP, kappa, TempE, error, TempF

N     = 50
TempF = 25.0d0
kappa = 0.3d0
Temp0 = 80.0d0
TempA = 20.0d0
dt    = TempF / dble(N)

t    = 0.0d0
Temp = Temp0

open(10,file="NLawC.dat")

do k = 0,N

    TempE = TempA + (Temp0 - TempA)*exp(-kappa*t)
    error = Temp - TempE

    print *, t, Temp
    write(10,*) t, Temp, TempE, error

    TempP = Temp
    Temp  = TempP + kappa*(TempA - TempP)*dt

    t = t + dt

end do

close(10)

end program 
