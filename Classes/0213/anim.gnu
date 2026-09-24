set terminal gif animate delay 5 
set output 'gaussian2D_time.gif'

set size square

set xlabel "{/Italic=25 x}"
set ylabel "{/Italic=25 y}"
set zlabel "{/Italic=25 f}"

set title '{/Italic=25 Gaussian Oscilation}'

do for[j=0:100]{

set contour
set xtics 1.0

splot [-5:5][-5:5][-1.1:1.1] \
'gaussian2D_time.dat' i j t '{/Italic=15 Scalar field}' w l

}
