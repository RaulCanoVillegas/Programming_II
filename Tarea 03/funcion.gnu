# Versión premium para paper con raíces destacadas
set terminal pdfcairo enhanced color size 7,5 font "Latin Modern Roman,12"
set output "funcion_paper_premium.pdf"

# Color personalizado
color_personal = "#286969"  # RGB {40,105,105}

# Configuración de estilo profesional
set style line 101 lc rgb '#000000' lt 1 lw 1.5  # Ejes
set style line 102 lc rgb color_personal lt 1 lw 3  # Función
set style line 103 lc rgb '#d95f02' lt 1 lw 2 pt 7 ps 1.5  # Raíces
set style line 104 lc rgb '#666666' lt 2 lw 1 dt 2  # Líneas auxiliares

# Configuración de bordes
set border 3 back ls 101
set tics nomirror
set grid back ls 104

# Títulos con formato LaTeX
set title "Gráfica para f(x) = (x-2)^4 - 10(x-2)^2" font "Latin Modern Roman,14"
set xlabel "x" font "Latin Modern Roman,12"
set ylabel "y" font "Latin Modern Roman,12"

# Rango y muestreo
set xrange [-2.5:6.5]
set yrange [-30:30]
set samples 2000

# Personalización de marcas
set xtics 1
set ytics 10
set mxtics 2
set mytics 2

# Calcular raíces
x1 = 2 - sqrt(10)  # ≈ -1.1623
x2 = 2              # = 2
x3 = 2 + sqrt(10)  # ≈ 5.1623

# Función
f(x) = (x-2)**4 - 10*(x-2)**2

# Crear puntos para las raíces
set print $roots
print sprintf("%.6f 0", x1)
print sprintf("%.6f 0", x2)
print sprintf("%.6f 0", x3)
set print

# Añadir anotaciones de las raíces con tamaño más grande (14 en lugar de 11)
set label 1 sprintf("x_3 = %.3f", x1) at x1-0.5, 2 font "Latin Modern Roman,14"  # Aumentado de 11 a 14
set label 2 sprintf("x_{1,2} = %.3f", x2) at x2+0.3, 2 font "Latin Modern Roman,14"  # Aumentado de 11 a 14
set label 3 sprintf("x_4 = %.3f", x3) at x3+0.3, 2 font "Latin Modern Roman,14"  # Aumentado de 11 a 14

# Graficar
plot f(x) w l ls 102 title "f(x)", \
     $roots using 1:2 w p ls 103 title "Raíces", \
     $roots using 1:2:(0):(-f($1)) w vectors heads size 0.1,20,60 ls 104 notitle

# Versión PNG
set terminal pngcairo enhanced size 1400,1000
set output "funcion_paper_premium.png"
replot

# Mostrar información de las raíces en consola
print "Raíces de la función:"
print sprintf("x₁ = %.6f", x1)
print sprintf("x₂ = %.6f", x2)
print sprintf("x₃ = %.6f", x3)
