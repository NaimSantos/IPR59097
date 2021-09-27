
set xrange[0:0.25]
plot "numeric_data_k_19_h_380.txt" using 1:2 title "k = 19" w l lw 2 lc "red", "numeric_data_k_59_h_380.txt" using 1:2 title "k = 59" w l lw 2 lc "blue", "numeric_data_k_210_h_380.txt" using 1:2 title "k = 210" w l lw 2 lc "yellow","numeric_data_k_400_h_380.txt" using 1:2 title "k = 400" w l lw 2 lc "forest-green",
set terminal pngcairo
set output "Grafico_h_fixo.png"
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