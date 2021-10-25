reset
set ter pdfcairo size 10 cm, 10 cm

set style line 1 lc 8 pt 7 
set style line 2 lc 2 pt 6 

set out "meas.pdf"

set yrange [0:1]

plot 'meas.txt' u 0:1 ls 1 w lines title 'plaquette', 'meas.txt' u 0:2 ls 2 w lines title 'acceptance rate'

set out
