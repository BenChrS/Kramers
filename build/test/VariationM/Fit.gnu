set terminal postscript 30 enhanced eps color dl 4
set output 'skalenverhalten.eps'

set xrange[0:4]
set yrange[0:10]

set xlabel 'm'
set ylabel 'R(m)'

f(x)=a*x**b
g(x)=c*x**d

c=0.1
d=-0.8

h=0.327877*0.135335*0.1*sqrt(2)

#fit [1/16:4] f(x) 'skalenverhalten.txt' using 1:($2/h) via a,b 
fit [1/16:4] g(x) 'skalenverhalten.txt' using 1:($3/h) via c,d

#plot \
#'skalenverhalten.txt' using 1:($2/h) title 'R_{k0}' w p lw 5, 'skalenverhalten.txt' using 1:($3/h) title 'R_D' w p lw 5,f(x) w l lw 2, g(x) w l lw 2

plot 'skalenverhalten.txt' using 1:($3/h) title 'R_D' w p lw 5,g(x) w l lw 2