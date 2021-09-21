L = 0.25
m2 = 6440.68
m = 80.2538
Tamb = 20
T0 = 320
TL = 75
f(x) = 300*(((55/300)*sinh(m*x)) + sinh(m*(L - x)))/ (sinh(m*L)) + 20

set xrange [0:0.25]
plot f(x) title 'Analítico' w l lw 2 lc "blue", "temperature_out1.dat" using 1:2 title "Numérico" lw 2 lc "red" dt '-'
set terminal pngcairo
set output "Graf3.png"
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