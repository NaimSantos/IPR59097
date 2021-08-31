set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 400, 400 
# set terminal pngcairo transparent enhanced font "arial,10" fontscale 1.0 size 400, 400 
set output 'heatmap.png'
unset key
set view map scale 1
set style data lines
set xtics border in scale 0,0 mirror norotate  autojustify
set ytics border in scale 0,0 mirror norotate  autojustify
set ztics border in scale 0,0 nomirror norotate  autojustify
unset cbtics
set rtics axis in scale 0,0 nomirror norotate  autojustify
set title "Perfil de Temperatura na placa" 
## set xrange [ -0.2500000 : 0.750000 ]
## set yrange [ -0.2500000 : 0.750000 ]

## set zrange [ * : * ] noreverse writeback
set cblabel "Escala de Temperatura" 
## set rrange [ * : * ] noreverse writeback
set palette rgbformulae 33, 13, 10
## NO_ANIMATION = 1
splot 'Dados.dat' matrix with image