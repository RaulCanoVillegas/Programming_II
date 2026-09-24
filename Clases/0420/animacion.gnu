# animacion.gnu

set terminal qt enhanced font "Arial,12"

set xlabel "x"
set ylabel "{/Symbol f}"
set xrange [-1:1]
set yrange [-1.5:1.5]
set grid
set key top right

do for [i=0:200] {
    set title sprintf("frame %d", i)
    plot 'wave_200.dat' index i using 1:3 with lines lw 2 lc rgb "blue" title "{/Symbol f}(x,t)"
    pause 0.05
}
