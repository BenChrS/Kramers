#plot
set terminal postscript 30 enhanced eps color dl 4
set output 'pot.eps'

set key b r 
set title 'Potenzial'
set xrange [0:3]
set yrange [0:4]
set xlabel 'V(x)'
set ylabel 'x'
plot \
"potential(0).txt" using 1:2 title 'pot' w l ls 2


