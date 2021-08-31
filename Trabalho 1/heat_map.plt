set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 400, 400
set output 'heatmap.png'
unset key
set view map scale 1
set xtics ("0" 1, "0.01" 10, "0.02" 20, "0.03" 30, "0.04" 40, "0.05" 49)
set ytics ("0" 1, "0.01" 10, "0.02" 20, "0.03" 30, "0.04" 40, "0.05" 49)
set xlabel "x (m)"
set ylabel "y (m)"
set xtics border in scale 0,0 mirror autojustify
set ytics border in scale 0,0 mirror autojustify
#set ztics border in scale 0,0 nomirror norotate  autojustify
# unset cbtics # remove a escala numerica de cd
#set rtics axis in scale 0,0 nomirror norotate  autojustify
set title "Perfil de Temperatura na placa" 
set cblabel "Temperatura" 
#set palette rgbformulae -34, -35, -36
set palette defined (100 "yellow", 101.5 "orange", 103 "tan1", 105 "red" )
splot 'Dados.dat' matrix with image