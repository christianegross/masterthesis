reset
set ter pdfcairo size 10 cm, 10 cm

set style line 1 lc 7 pt 7 
set style line 2 lc 2 pt 6 

file='measbeta1.0Nt16Ns16cnum.txt'
file2='measbeta1.0Nt16Ns16.txt'

set out 'averageplaquette.pdf'

set title 'Delta tau=1, N_s=N_t=8'
set ylabel '<P>'
set xlabel 'beta'
plot 'meascombined.txt' u 1:5 ls 1 title ''

set out "meas.pdf"

set yrange [0:1]

set xlabel 'configuration'
set ylabel ''

plot file u 0:1 ls 1 w lines title 'plaquette', file u 0:2 ls 2 w lines title 'acceptance rate'
plot file2 u 0:1 ls 1 w lines title 'plaquette', file2 u 0:2 ls 2 w lines title 'acceptance rate'


set yrange [0.645:0.65]
plot file u 0:1 ls 2 ps 0.3 title '', file2 u 0:1 ls 1 ps 0.1 title ''

set out
