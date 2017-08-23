#plot
set terminal postscript 30 enhanced eps color dl 4
set output 'Ekin.eps'

set key b c 
set title 'kinetisches Energie'
#set xrange [-2:2]
#set yrange [0:180]
set xlabel 't'
set ylabel 'E_{kin}(t)'
plot \
"kinEnergyAv(0).txt" using 1:2 title 'numerisch' w l ls 2, "kinEnergyTheo(0).txt" using 1:2 title 'analytisch' w l ls 1


