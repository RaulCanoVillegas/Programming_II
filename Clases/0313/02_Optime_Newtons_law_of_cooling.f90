program Newtons_law_of_cooling
implicit none

integer :: k,N
real(kind=8) :: t, dt, TempA, Temp0, Temp, TempP, kappa, TempE, error, TempF, k1, k2

N     = 50
TempF = 25.0d0
kappa = 0.3d0
Temp0 = 80.0d0
TempA = 20.0d0
dt    = TempF / dble(N)

t     = 0.0d0
Temp  = Temp0
TempE = Temp0
error = Temp - TempE

open(10,file="NLawC_Optime.dat")

do k = 0,N

    print *, t, Temp, TempE, error
    TempP = Temp
    !--------------------------------------------
    k1    = kappa * (TempA - TempP)
    Temp  = k1 * dt + TempP
    k2    = kappa * (TempA - Temp ) 
    Temp  = k2 * dt - TempP
    !--------------------------------------------
    t     = t + dt
    TempE = TempA + (Temp0 - TempA)*exp(-kappa*t)

    error = Temp - TempE

end do

close(10)

end program 
