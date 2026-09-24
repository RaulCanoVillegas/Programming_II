program rk_newtons_law_of_cooling
implicit none

integer :: k,N
real(kind=8) :: t, dt, TempA, Temp0, Temp, TempP
real(kind=8) :: kappa, TempE, error, TempF, k1, k2

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

open(10,file="NLawC_RK.dat")

do k = 0,N

    TempE = TempA + (Temp0 - TempA)*exp(-kappa*t)
    error = Temp - TempE

    print *, t, Temp, TempE, error
    write(10,*) t, Temp, TempE, error

    TempP = Temp

    ! Runge-Kutta orden 2
    k1 = kappa*(TempA - TempP)
    k2 = kappa*(TempA - (TempP + dt*k1))
    Temp = TempP + 0.5d0*(k1 + k2)*dt

    t = t + dt

end do

close(10)

end program
