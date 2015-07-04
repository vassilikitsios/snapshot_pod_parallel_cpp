#
# gnuplot> load 'plot_snapshot_pod_results_eps.gp'
#


# -------------------------------- Ouput settings -------------------------------- 
set terminal postscript enhanced font "Times-Roman" 20
set output 'view_output.eps'


# -------------------------------- Eigenvalues -------------------------------- 

set title 'Snapshot POD eigenvalues'
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

set title 'Snapshot POD temporal modes 1 and 2'
plot '../results/temporal_mode_0001.dat' using 1:2 title "mode 1" with lines lt 1,	'../results/temporal_mode_0002.dat' using 1:2 title "mode 2" with lines lt 2

set title 'Snapshot POD temporal modes 3 and 4'
plot '../results/temporal_mode_0003.dat' using 1:2 title "mode 3" with lines lt 1,	'../results/temporal_mode_0004.dat' using 1:2 title "mode 4" with lines lt 2

set title 'Snapshot POD temporal modes 5 and 6'
plot '../results/temporal_mode_0005.dat' using 1:2 title "mode 5" with lines lt 1,	'../results/temporal_mode_0006.dat' using 1:2 title "mode 6" with lines lt 2

set title 'Snapshot POD temporal modes 7 and 8'
plot '../results/temporal_mode_0007.dat' using 1:2 title "mode 7" with lines lt 1,	'../results/temporal_mode_0008.dat' using 1:2 title "mode 8" with lines lt 2

set title 'Snapshot POD temporal modes 9 and 10'
plot '../results/temporal_mode_0009.dat' using 1:2 title "mode 9" with lines lt 1,	'../results/temporal_mode_0010.dat' using 1:2 title "mode 10" with lines lt 2

# -------------------------------- EOF -------------------------------- 
