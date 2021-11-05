reset
set ter pdfcairo size 10 cm, 10 cm
set samples 2000
set style line 1 lc 7 pt 7 
set style line 2 lc 2 pt 6 
set style line 3 lc 3 pt 5
set style line 4 lc 4 pt 4
set style line 5 lc 5 pt 3

filesu2='meas/output.u1-metropolis.data'

d=4
weakcoupling(x)=1.0-1.0/(4.0*x)-1.0/(32.0*x**2)-0.0131/x**3-0.00752/x**4
F(x)=d*(d-1)*(x**2/8.0)#-x**4/128.0+(d/192.0-11.0/1152.0)*x**6+(757.0/98304.0-d/256.0)*x**8*(d*d/512.0-85.0/12288.0*d+2473.0/409600.0)*x**10+(2467.0*d/262144.0-29*d*d/12288.0-1992533.0/212336640.0)*x**12+(5.0*d**3/4096.0-237.0*d**2/32768.0+178003.0*d/11796480.0-38197099.0/3468165120.0)*x**14+(-15/8192.0*d**3+1485.0/131072.0*d*d-53956913.0*d/2264924160.0+11483169709.0/676457349120.0)*x**16)
Pstrong(x)=2*(x/4-4*x**3/128.0+(d/192.0-11.0/1152.0)*6*x**5+(757.0/98304.0-d/256.0)*8*x**7*(d*d/512.0-85.0/12288.0*d+2473.0/409600.0)*10*x**9+(2467.0*d/262144.0-29*d*d/12288.0-1992533.0/212336640.0)*12*x**11+(5.0*d**3/4096.0-237.0*d**2/32768.0+178003.0*d/11796480.0-38197099.0/3468165120.0)*14*x**13+(-15/8192.0*d**3+1485.0/131072.0*d*d-53956913.0*d/2264924160.0+11483169709.0/676457349120.0)*16*x**15)
#*d*(d-1)

set out 'averageplaquette.pdf'
binsize=72

set title 'Delta tau=1, N_s=N_t=10, binsize='.binsize
set ylabel '<1-P>'
set xlabel 'beta'
set key bottom right

set yrange [0:1]
set xrange[0:6]
plot 'measplaquette.txt' u ($5==0?$1:1/0):(1-$6):7 w yerrorbars ls 1 ps 0.2 title 'naive error',\
'measplaquette.txt' u ($5==binsize?$1:1/0):(1-$6):7 w yerrorbars ls 4 ps 0.2 title 'errors C',\
 weakcoupling(x) ls 2 title 'weak', Pstrong(x) ls 3 title '2*dF/dbeta*1/(d*(d-1))',\
 'meas/resulthadron.csv' u ($2==binsize?$8:1/0):(1-$3):4 w yerrorbars ls 5 ps 0.1 title 'errors hadron'
set xrange[0:1]
unset yrange
plot 'measplaquette.txt' u ($5==0?$1:1/0):(1-($6)-($1/2.0)):7 w yerrorbars ls 1 ps 0.2 title 'naive error',\
'measplaquette.txt' u ($5==binsize?$1:1/0):(1-($6)-($1/2.0)):7 w yerrorbars ls 4 ps 0.2 title 'errors C',\
 Pstrong(x)-x/2.0 ls 3 title '2*dF/dbeta*1/(d*(d-1))-beta/2',\
 'meas/resulthadron.csv' u ($2==binsize?$8:1/0):(1-$3-($8/2)):4 w yerrorbars ls 5 ps 0.1 title 'errors hadron'
unset yrange
unset xrange

set ylabel "error <P>"
set key c tm
#do for [elem in "2, 4, 8, 16, 32, 52, 72, 92, 112, 132, 152, 172, 192, 212, 232"] {
#len=int(elem)
#set title "len=".elem
#plot 'measplaquette.txt' u (($5==len)?$1:1/0):7 ls 1 title 'C', 'meas/resulthadron.csv' u (($2==len)?$8:1/0):4:5 w yerrorbars ls 2 title 'hadron'
#}

set xrange[0.5:233]
do for [beta in "0.05 0.1 0.2 0.3 0.4 0.6 0.8 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 1.2 1.4 1.6 1.8"]{
set title "beta=".beta
plot 'measplaquette.txt' u (($1==beta)?$5:1/0):7 ls 1 title 'C', 'meas/resulthadron.csv' u (($8==beta)?$2:1/0):4:5 w yerrorbars ls 2 title 'hadron'
}

 
#set xrange[3:100]
#unset yrange
##plot 'meascombined.txt' u 1:(1-$5):6 w yerrorbars ls 1 ps 0.2 title 'measured', weakcoupling(x) ls 2 title 'weak'#, F(x) ls 3 title 'strong'
#unset xrange 
#set logscale x
#set xlabel 'binsize'
#set ylabel 'std(P)'
#plot 'measplaquette.txt' u 5:7 ls 1 title ''
#unset logscale x

set title 'beta=sqrt(2), N_s=N_t=8'
#set ylabel '<P>'
set xrange [0.05:1.05]
set yrange [0.1:0.3]
set xlabel 'Delta tau'
#plot 'plaquettedeltatau.txt' u 2:5:6 w yerrorbars ls 1 ps 0.2 title 'plaquette'

unset title

set out "comparison.pdf"
unset xrange
set yrange [0.8:1]
set xlabel 'step'
set ylabel '<Re(Tr(U_{mu}(x)U_{nu}(x+mu)U^{dagger}_{mu}(x+nu)U^{dagger}_{nu}(x)>'

#plot filesu2 u 1:2 ls 1 w lines title 'SU2-code', 'meas/measbeta2.100000Nt10Ns10meas1000.txt' u 0:(1-$1) ls 2 w lines title 'U1-code'
set xlabel 'step'
set ylabel '<1-P>'

set out "meas.pdf"
#plot filesu2 u 1:2 ls 1 w lines title 'SU2-code', 'meas/measbeta2.100000Nt10Ns10meas1000.txt' u 0:(1-$1) ls 2 w lines title 'U1-code'

set yrange [0:1]
plot 'meas/thermbeta5.500000Nt10Ns10.txt' u 0:(1-$1) ls 1 w lines title 'beta=5.5', 'meas/measbeta5.500000Nt10Ns10.txt' u ($0+5000):(1-$1) ls 2 w lines title 'beta=5.5'
plot 'meas/thermbeta0.400000Nt10Ns10.txt' u 0:(1-$1) ls 1 w lines title 'beta=0.4', 'meas/measbeta0.400000Nt10Ns10.txt' u ($0+5000):(1-$1) ls 2 w lines title 'beta=0.4'
plot 'meas/thermbeta0.050000Nt10Ns10.txt' u 0:(1-$1) ls 1 w lines title 'beta=0.05', 'meas/measbeta0.050000Nt10Ns10.txt' u ($0+5000):(1-$1) ls 2 w lines title 'beta=0.05'

set out
