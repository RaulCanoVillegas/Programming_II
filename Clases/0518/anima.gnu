# ============================================================
#  anim3D.gnu  –  Animacion 3D (superficie) de la onda 2D
#  Uso:  gnuplot anim3D.gnu
# ============================================================

Nt_total = 800
cada     = 10
Nframes  = int(Nt_total / cada)

set terminal gif animate delay 5 size 800,700 optimize
set output 'wave2D_3D.gif'

set palette defined (0 "navy", 0.3 "cyan", 0.5 "white", 0.7 "yellow", 1 "red")
set cbrange [-0.5:1]

set xlabel "x" font ",10"
set ylabel "y" font ",10"
set zlabel "phi" font ",10"
set xrange [-1:1]
set yrange [-1:1]
set zrange [-0.5:1]
set hidden3d
unset key

# Angulo de vista
set view 55, 30

do for [k=0:Nframes] {
    set title sprintf("Onda 2D – superficie  (frame %03d)", k) font ",12"

    splot 'EcOnda2D.dat' index k \
          using 1:2:3 with pm3d \
          title ''
}

set output
print "GIF 3D generado: onda2D_3D.gif"
