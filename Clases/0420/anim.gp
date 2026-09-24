set yrange[-1:1]
set xrange[-1:1]

do for [i=0:*] {
    plot 'wave.dat' index i using 1:3 with lines title sprintf("bloque %d", i)
}
