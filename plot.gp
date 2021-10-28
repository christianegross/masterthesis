reset
set ter pdfcairo size 10 cm, 10 cm

set style line 1 lc 7 pt 7 
set style line 2 lc 2 pt 6 
set style line 3 lc 3 pt 5

file='measbeta1.0Nt16Ns16cnum.txt'
file2='measbeta1.0Nt16Ns16.txt'
file3='measbeta5.0Nt8Ns8random.txt'

d=4
weakcoupling(x)=1.0-1.0/(4.0*x)-1.0/(32.0*x**2)-0.0131/x**3-0.00752/x**4
F(x)=d*(d-1)*(x**2/8.0-x**4/128.0+(d/192.0-11.0/1152.0)*x**6+(757.0/98304.0-d/256.0)*x**8*(d*d/512.0-85.0/12288.0*d+2473.0/409600.0)*x**10+(2467.0*d/262144.0-29*d*d/12288.0-1992533.0/212336640.0)*x**12+(5.0*d**3/4096.0-237.0*d**2/32768.0+178003.0*d/11796480.0-38197099.0/3468165120.0)*x**14+(-15/8192.0*d**3+1485.0/131072.0*d*d-53956913.0*d/2264924160.0+11483169709.0/676457349120.0)*x**16)
#*d*(d-1)

set out 'averageplaquette.pdf'

set title 'Delta tau=1, N_s=N_t=8'
set ylabel '<P>'
set xlabel 'beta'
set key bottom right

set yrange [0:1]
set xrange[0.05:3]
plot 'meascombined.txt' u 1:(1-$5):6 w yerrorbars ls 1 ps 0.2 title 'measured', weakcoupling(x) ls 2 title 'weak'#, F(x) ls 3 title 'strong'
set xrange[3:100]
unset yrange
plot 'meascombined.txt' u 1:(1-$5):6 w yerrorbars ls 1 ps 0.2 title 'measured', weakcoupling(x) ls 2 title 'weak'#, F(x) ls 3 title 'strong'
unset title

set out "meas.pdf"

set yrange [0.7:0.75]

plot 'measbeta1.4Nt8Ns8ordered.txt' u 0:1 ls 1 w lines title 'ordered', 'measbeta1.4Nt8Ns8random.txt' u 0:1 ls 2 w lines title 'random'

set yrange [0:1]

set xlabel 'configuration'
set ylabel ''

plot file u 0:1 ls 1 w lines title 'plaquette', file u 0:2 ls 2 w lines title 'acceptance rate'
plot file2 u 0:1 ls 1 w lines title 'plaquette', file2 u 0:2 ls 2 w lines title 'acceptance rate'
plot file3 u 0:1 ls 1 w lines title 'plaquette', file3 u 0:2 ls 2 w lines title 'acceptance rate'


set yrange [0.645:0.65]
plot file u 0:1 ls 2 ps 0.3 title '', file2 u 0:1 ls 1 ps 0.1 title ''

set out
