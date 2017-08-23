#plot
set terminal postscript 30 enhanced eps color dl 4
set output 'fluss.eps'

set key b r 
set title 'Fluss'
set xrange [0:100]
#set yrange [0:1]
set xlabel 't'
set ylabel 'j_+(t)'
f(x)=0.00463191
plot \
"fluxPositiveTotal(0).txt" using ($1*6.44636):2 title 'n=800000,dt=0.012207' w l ls 2,f(x) w l ls 3




