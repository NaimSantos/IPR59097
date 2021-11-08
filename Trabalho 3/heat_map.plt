set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 500, 500
set output 'heatmap.png'
unset key
set view map scale 1
set xtics ("0" 1, "0.01" 20, "0.02" 40, "0.03" 60, "0.04" 80, "0.05" 99)
set ytics ("0" 1, "0.01" 20, "0.02" 40, "0.03" 60, "0.04" 80, "0.05" 99)
set xlabel "x (m)"
set ylabel "y (m)"
set xtics border in scale 0,0 mirror autojustify
set ytics border in scale 0,0 mirror autojustify
#set ztics border in scale 0,0 nomirror norotate  autojustify
# unset cbtics # remove a escala numerica de cd
#set rtics axis in scale 0,0 nomirror norotate  autojustify
set title "Perfil de Temperatura na placa" 
set cblabel "Temperatura" 
# set palette rgbformulae -34, -35, -36
set palette defined (20 'white', 21 "yellow", 40 "orange", 60 "tan1", 90 "red", 120 "dark-red" )
#set palette defined (100 "yellow", 101.5 "orange", 103 "tan1", 105 "red" )
splot '2D_Heat.txt' matrix with image