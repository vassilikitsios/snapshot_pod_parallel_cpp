
# -------------------------------- Ouput settings -------------------------------- 
set terminal postscript enhanced eps font "Times-Roman" 32 dashlength 3

# -------------------------------- Eigenvalues -------------------------------- 

set output './eigenvalues.eps'
set logscale y 10
set logscale y2 10
set xlabel 'mode'
set ylabel 'energy' rotate by 90
set y2label 'cumulative energy' rotate by 90
set ytics nomirror
set y2tics nomirror
plot '../results/eigenvalues.dat' using 1:3 title "eigenvalues" with linespoints lt 1 pt 1 axes x1y1, \
                  '../results/eigenvalues.dat' using 1:4 title "cumulative" with linespoints lt 2 pt 2 axes x1y2

# -------------------------------- Temporal Modes -------------------------------- 

unset logscale
set xlabel 'time'
set ylabel 'fourier coefficients' rotate by 90
unset y2label
unset y2tics
set ytics mirror

set output './temporal_modes_0001_0002.eps'
plot '../results/temporal_mode_0001.dat' using 1:2 title "mode 1" with lines lt 1,	'../results/temporal_mode_0002.dat' using 1:2 title "mode 2" with lines lt 2

set output './temporal_modes_0003_0004.eps'
plot '../results/temporal_mode_0003.dat' using 1:2 title "mode 3" with lines lt 1,	'../results/temporal_mode_0004.dat' using 1:2 title "mode 4" with lines lt 2

set output './temporal_modes_0005_0006.eps'
plot '../results/temporal_mode_0005.dat' using 1:2 title "mode 5" with lines lt 1,	'../results/temporal_mode_0006.dat' using 1:2 title "mode 6" with lines lt 2

set output './temporal_modes_0007_0008.eps'
plot '../results/temporal_mode_0007.dat' using 1:2 title "mode 7" with lines lt 1,	'../results/temporal_mode_0008.dat' using 1:2 title "mode 8" with lines lt 2

set output './temporal_modes_0009_0010.eps'
plot '../results/temporal_mode_0009.dat' using 1:2 title "mode 9" with lines lt 1,	'../results/temporal_mode_0010.dat' using 1:2 title "mode 10" with lines lt 2

# -------------------------------- EOF -------------------------------- 
