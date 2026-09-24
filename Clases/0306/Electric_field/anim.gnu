set terminal gif animate delay 5
set output 'electric_field.gif'

set size square

set xlabel 'x'
set ylabel 'y'

set xrange[-1:1]
set yrange[-1:1]

do for [i=0:100] {

    plot 'animated_field.dat' index i \
    u 1:2:3:4 w vectors head filled lw 1 \
    t 'Electric field'

}
