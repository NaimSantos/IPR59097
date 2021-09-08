plot "temperature_out.txt" using 1:2 notitle w l lw 2 lc "red"
set terminal pngcairo
set output "Graf1.png"
set grid
set title "Perfil de Temperatura na aleta"
set xlabel "Posição (m)"
set ylabel "Temperatura ({\260}C)"
replot
set terminal wxt
set output


# with lines = w l
# linewidth = lw
# linecolor = lc
# dashtype = dt