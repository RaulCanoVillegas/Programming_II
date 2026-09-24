# ============================================================
#  anim3D.gnu  –  Animacion 3D simple (malla) de la onda 2D
#  Uso:  gnuplot anim3D.gnu
# ============================================================

Nt_total = 800
cada     = 10
Nframes  = int(Nt_total / cada)

set terminal gif animate delay 5 size 700,600 optimize
set output 'onda2D_3D.gif'

set xlabel "x"
set ylabel "y"
set zlabel "phi"
set xrange [-1:1]
set yrange [-1:1]
set zrange [-0.5:1]
set hidden3d
unset key
unset colorbox

# Vista similar a la imagen de referencia
set view 60, 30

do for [k=0:Nframes] {
    set title sprintf("Ecuacion de onda 2D  (frame %03d)", k)

    splot 'EcOnda2D.dat' index k \
          using 1:2:3 with lines \
          lc rgb "purple" lw 0.5
}

set output
print "GIF generado: onda2D.gif"
