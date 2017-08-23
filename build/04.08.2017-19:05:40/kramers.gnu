#plot
set terminal postscript 30 enhanced eps color dl 4
set output 'fluss.eps'

set key t r 
set title 'Fluss'
set xrange [0:80]
set yrange [0:0.0001]
set xlabel 't'
set ylabel 'j_+(t)'
k(x)=0.0000210514
plot \
"fluxPositiveTotal(0).txt" using ($1*6.44636):2 title 'n=80000,dt=0.0195312' w l ls 2,k(x) w l ls 3




