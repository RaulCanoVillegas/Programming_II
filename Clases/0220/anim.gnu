set terminal gif animate delay 5 
set output 'animacion.gif'

set size square

set xlabel "{/Italic=25 x}"
set ylabel "{/Italic=25 y}"
set title '{/Italic=25 Vector field}'

unset key
do for[j=0:100]{

#set contour
#set xtics 1.0

plot [-1.5:1.5][-1.5:1.5] 'vector_field.dat' u 1:2:3:4 i j w vectors

}

