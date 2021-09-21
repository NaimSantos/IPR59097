
set xrange[0:0.25]
plot "analitico.txt" using 1:2 title 'Analítico' w l lw 2 lc "blue", "dataN51.txt" using 1:2 title "Numérico" lw 2 lc "red" dt '-'
set terminal pngcairo
set output "Graf5.png"
set grid
#set title "Perfil de Temperatura na aleta"
set xlabel "Posição (m)"
set ylabel "Temperatura ({\260}C)"
replot
set terminal wxt
set output


# with lines = w l
# linewidth = lw
# linecolor = lc
# dashtype = dt