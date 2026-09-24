# ============================================================
#  animate_wave.gp  —  Animacion de la ecuacion de onda 1D
#  Uso:  gnuplot animate_wave.gp
# ============================================================

datafile = "wave.dat"

# --- Parametros de la ventana ---
set terminal qt size 900,500 title "Ecuacion de onda 1D — RK2"
set xlabel "x"
set ylabel "{/Symbol f}(x,t)"
set xrange [-1:1]
set yrange [-1.2:1.2]
set grid lc rgb "#cccccc"
set key off
set style line 1 lc rgb "#1a7abf" lw 2

# --- Contar cuantos bloques (frames) hay en el archivo ---
#     Cada bloque esta separado por dos lineas en blanco.
stats datafile nooutput
nframes = int(STATS_blocks)

print sprintf("Frames encontrados: %d", nframes)

# --- Tiempo entre frames (ajusta a tu gusto) ---
pause_time = 0.04          # segundos entre frames

# --- Animacion: un loop sobre los bloques del archivo ---
do for [i=0:nframes-1] {
    set title sprintf("Frame %d / %d", i, nframes-1) font ",12"
    plot datafile index i \
         using 1:2 with lines ls 1 title ""
    pause pause_time
}

# --- Al terminar, mostrar el ultimo frame y esperar ---
set title sprintf("Frame final (%d)", nframes-1) font ",12"
plot datafile index (nframes-1) using 1:2 with lines ls 1
pause -1 "Presiona Enter para salir..."
