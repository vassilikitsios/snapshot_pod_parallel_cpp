
# -------------------------------- Ouput settings -------------------------------- 
set terminal postscript enhanced eps font "Times-Roman" 32 dashlength 3

# -------------------------------- Max Reconstriction L2 Error -------------------------------- 

set logscale y 10
set xlabel 'mode'
set ylabel 'max L2 Error' rotate by 90
set output './max_spatial_reconstruction_L2_error.eps'
plot '../results/max_spatial_reconstruction_L2_error.dat' using 1:2 title "u" with linespoints lt 2 pt 2, \
	'../results/max_spatial_reconstruction_L2_error.dat' using 1:3 title "v" with linespoints lt 3 pt 3, \
	'../results/max_spatial_reconstruction_L2_error.dat' using 1:4 title "w" with linespoints lt 4 pt 4, \
	'../results/max_spatial_reconstruction_L2_error.dat' using 1:5 title "pressure" with linespoints lt 5 pt 5, \
	'../results/max_spatial_reconstruction_L2_error.dat' using 1:6 title "total" with linespoints lt 6 pt 6


# -------------------------------- Reconstriction L2 Error -------------------------------- 

unset logscale
set xlabel 'time'
set ylabel 'L2 Error' rotate by 90

set output './spatial_reconstruction_L2_error_0001.eps'
plot '../results/spatial_reconstruction_L2_error_0001.dat' using 1:2 title "u" with lines lt 2,	\
	'../results/spatial_reconstruction_L2_error_0001.dat' using 1:3 title "v" with lines lt 3, \
	'../results/spatial_reconstruction_L2_error_0001.dat' using 1:4 title "w" with lines lt 4, \
	'../results/spatial_reconstruction_L2_error_0001.dat' using 1:5 title "p" with lines lt 5, \
	'../results/spatial_reconstruction_L2_error_0001.dat' using 1:6 title "total" with lines lt 6

set output './spatial_reconstruction_L2_error_0002.eps'
plot '../results/spatial_reconstruction_L2_error_0002.dat' using 1:2 title "u" with lines lt 2,	\
	'../results/spatial_reconstruction_L2_error_0002.dat' using 1:3 title "v" with lines lt 3, \
	'../results/spatial_reconstruction_L2_error_0002.dat' using 1:4 title "w" with lines lt 4, \
	'../results/spatial_reconstruction_L2_error_0002.dat' using 1:5 title "p" with lines lt 5, \
	'../results/spatial_reconstruction_L2_error_0002.dat' using 1:6 title "total" with lines lt 6



# -------------------------------- EOF -------------------------------- 
