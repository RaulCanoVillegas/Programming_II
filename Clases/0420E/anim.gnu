set terminal gif animate delay 5
set output 'anim.gif'

#set size square
set xlabel "{/Italic=25 x}"
set ylabel "{/Italic=25 u}"

#set tittle '{/Italic=25 Campo Vectorial}'
unset key
do for[j=0:400]{

#set contour
#set xtics 1.0
plot [-1:1][-1.1:1.1] 'Wave' i j w l lw 3

}
