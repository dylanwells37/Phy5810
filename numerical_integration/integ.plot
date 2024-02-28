# gnuplot file

set xlabel "log10(N)"

set ylabel "log10(relative error)"

set title "Relative error vs N"

# fit a line to milne and simpson from -2 to -0.5
f(x) = a*x + b
f2(x) = a2*x + b2
# range for fit -2 to -0.5
fit [2:3.5] f(x) "milne.dat" using 1:3 via a,b
fit [2:4] f2(x) "simpson.dat" using 1:3 via a2,b2

f3(x) = a3*x + b3
fit [0:2] f3(x) "milne.dat" using 1:3 via a3,b3
f4(x) = a4*x + b4
fit [0.0:2] f4(x) "simpson.dat" using 1:3 via a4,b4

f5(x) = a5*x + b5
fit [4:8.5] f5(x) "milne.dat" using 1:3 via a5,b5

fit_title = sprintf("Milne, %-+4.1f*x %-+4.1f",a,b)
fit_title2 = sprintf("Simpson, %-+4.1f*x %-+4.1f",a2,b2)
fit_title3 = sprintf("Milne, %-+4.1f*x %-+4.1f",a3,b3)
fit_title4 = sprintf("Simpson, %-+4.1f*x %-+4.1f",a4,b4)
fit_title5 = sprintf("Milne, %-+4.1f*x %-+4.1f",a5,b5)




plot "milne.dat" using 1:3 with lines title "Milne", \
    "simpson.dat" using 1:3 with lines title "Simpson", \
    "gsl.dat" using 1:3 with lines title "GSL", \
    f(x) title fit_title, \
    f2(x) title fit_title2, \
    f3(x) title fit_title3, \
    f4(x) title fit_title4, \
    f5(x) title fit_title5

# save to postscript
set term postscript eps color
set output "error.ps"
replot
