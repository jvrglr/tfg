set title "n(t) using logscale"
set key left
set grid
set logscale
set xlabel "t"
set ylabel "n(t)"
f(x)=b*(x**a)
g(x)=x**(0.313)
#g(x)=x**(-0.159)
#g(x)=x**(1.265)
fit f(x) "data_n.dat" u 1:2 via a,b
set term postscript eps enhanced color 14
set output "n(t)_11_2_17.eps"
plot "data_n.dat" u 1:2 w points pt 7 title "Simulated behavior",\
 f(x) title "Fit",g(x) title "Theoretical behavior (Marro)"
