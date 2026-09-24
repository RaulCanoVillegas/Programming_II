# Versión detallada para raíz triple
set terminal pdfcairo enhanced color size 8,6 font "Latin Modern Roman,12"
set output "funcion_cubica_detallada.pdf"

# Color personalizado
color_personal = "#286969"

# Estilos
set style line 101 lc rgb '#000000' lt 1 lw 1.5
set style line 102 lc rgb color_personal lt 1 lw 3
set style line 103 lc rgb '#d95f02' lt 1 lw 2 pt 7 ps 2.5
set style line 104 lc rgb '#666666' lt 2 lw 1 dt 2
set style line 105 lc rgb '#CC6600' lt 1 lw 2 dt 2
set style line 106 lc rgb '#0066CC' lt 1 lw 2 dt 3

set border 3 back ls 101
set tics nomirror
set grid back ls 104

set title "Análisis de raíz triple: f(x) = (x-2)^3" font "Latin Modern Roman,16"
set xlabel "x" font "Latin Modern Roman,14"
set ylabel "y" font "Latin Modern Roman,14"

set xrange [1:3]
set yrange [-1:1]
set samples 2000

set xtics 0.2
set ytics 0.2
set mxtics 4
set mytics 4

# Función
f(x) = (x-2)**3
df(x) = 3*(x-2)**2
ddf(x) = 6*(x-2)

x_raiz = 2.0

# Crear datos para la raíz
set print $raiz
print sprintf("%.6f 0", x_raiz)
set print

# Crear datos para la tangente
set print $tangente
do for [i=0:100] {
    x = 1.5 + i*0.01
    print sprintf("%.6f %.6f", x, df(x_raiz)*(x-x_raiz))
}
set print

# Anotaciones
set label 1 "Raíz triple" at x_raiz+0.15, 0.15 font "Latin Modern Roman,12"
set label 2 "f(x) = (x-2)^3" at 2.5, 0.8 font "Latin Modern Roman,12"
set label 3 "Tangente horizontal (f' = 0)" at 2.3, -0.3 font "Latin Modern Roman,10"
set label 4 "Punto de inflexión" at 1.7, -0.2 font "Latin Modern Roman,10"

# Dibujar línea vertical en la raíz
set arrow from x_raiz, -1 to x_raiz, 1 nohead ls 104

plot f(x) w l ls 102 title "f(x) = (x-2)^3", \
     $raiz using 1:2 w p ls 103 title "Raíz en x=2", \
     $tangente using 1:2 w l ls 105 title "Tangente (f'=0)"

# Versión PNG
set terminal pngcairo enhanced size 1400,1000
set output "funcion_cubica_detallada.png"
replot

print "Propiedades de la raíz triple en x=2:"
print sprintf("f(2) = %.6f", f(x_raiz))
print sprintf("f'(2) = %.6f", df(x_raiz))
print sprintf("f''(2) = %.6f", ddf(x_raiz))
print "La función tiene un punto de inflexión en la raíz"
