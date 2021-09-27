
set xrange[0:0.25]
plot "numeric_data_k_59_h_100.txt" using 1:2 title "h = 100" w l lw 2 lc "forest-green", "numeric_data_k_59_h_380.txt" using 1:2 title "h = 380" w l lw 2 lc "blue", "numeric_data_k_59_h_600.txt" using 1:2 title "h = 600" w l lw 2 lc "red"
set terminal pngcairo
set output "Grafico_k_fixo.png"
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