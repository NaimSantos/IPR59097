set terminal pngcairo transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
set output 'heatmaps.png'
unset key
set view map scale 1
set style data lines
set xtics border in scale 0,0 mirror norotate  autojustify
set ytics border in scale 0,0 mirror norotate  autojustify
set ztics border in scale 0,0 nomirror norotate  autojustify
unset cbtics
set rtics axis in scale 0,0 nomirror norotate  autojustify
set title "Perfil de Temperatura na placa" 
set xrange [ -0.2500000 : 0.750000 ] noreverse nowriteback
set x2range [ * : * ] noreverse writeback
set yrange [ -0.2500000 : 0.750000 ] noreverse nowriteback
set y2range [ * : * ] noreverse writeback
set zrange [ * : * ] noreverse writeback
set cblabel "Gradiente" 
set cbrange [ 0.00000 : 5.00000 ] noreverse nowriteback
set rrange [ * : * ] noreverse writeback
set palette rgbformulae -7, 2, -7
NO_ANIMATION = 1
## Last datafile plotted: "$map1"
splot 'Dados.dat' matrix with image