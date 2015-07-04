
# -------------------------------- Ouput settings -------------------------------- 
set terminal postscript enhanced font "Times-Roman" 20
set output 'view_output.eps'

# -------------------------------- L infinite Error -------------------------------- 

set xlabel 'time'

set title 'L_{\infty} Error for each Mode level reconstruction across all snapshots'
set logscale y 10
set xlabel 'mode'
set ylabel 'max L2 Error' rotate by 90
plot '../results/max_spatial_reconstruction_L2_error.dat' using 1:2 title "u" with linespoints lt 2 pt 2, \
	'../results/max_spatial_reconstruction_L2_error.dat' using 1:3 title "v" with linespoints lt 3 pt 3, \
	'../results/max_spatial_reconstruction_L2_error.dat' using 1:4 title "w" with linespoints lt 4 pt 4, \
	'../results/max_spatial_reconstruction_L2_error.dat' using 1:5 title "pressure" with linespoints lt 5 pt 5, \
	'../results/max_spatial_reconstruction_L2_error.dat' using 1:6 title "total" with linespoints lt 6 pt 6


# -------------------------------- Reconstriction L2 Error -------------------------------- 

unset logscale
set xlabel 'time'
set ylabel 'L2 Error' rotate by 90

set title 'Spatial Reconstruction L2 Error Using 1 Mode'
plot '../results/spatial_reconstruction_L2_error_0001.dat' using 1:2 title "u" with lines lt 2,	\
	'../results/spatial_reconstruction_L2_error_0001.dat' using 1:3 title "v" with lines lt 3, \
	'../results/spatial_reconstruction_L2_error_0001.dat' using 1:4 title "w" with lines lt 4, \
	'../results/spatial_reconstruction_L2_error_0001.dat' using 1:5 title "p" with lines lt 5, \
	'../results/spatial_reconstruction_L2_error_0001.dat' using 1:6 title "total" with lines lt 6

set title 'Spatial Reconstruction L2 Error Using 2 Modes'
plot '../results/spatial_reconstruction_L2_error_0002.dat' using 1:2 title "u" with lines lt 2,	\
	'../results/spatial_reconstruction_L2_error_0002.dat' using 1:3 title "v" with lines lt 3, \
	'../results/spatial_reconstruction_L2_error_0002.dat' using 1:4 title "w" with lines lt 4, \
	'../results/spatial_reconstruction_L2_error_0002.dat' using 1:5 title "p" with lines lt 5, \
	'../results/spatial_reconstruction_L2_error_0002.dat' using 1:6 title "total" with lines lt 6

set title 'Spatial Reconstruction L2 Error Using 3 Modes'
plot '../results/spatial_reconstruction_L2_error_0003.dat' using 1:2 title "u" with lines lt 2,	\
	'../results/spatial_reconstruction_L2_error_0003.dat' using 1:3 title "v" with lines lt 3, \
	'../results/spatial_reconstruction_L2_error_0003.dat' using 1:4 title "w" with lines lt 4, \
	'../results/spatial_reconstruction_L2_error_0003.dat' using 1:5 title "p" with lines lt 5, \
	'../results/spatial_reconstruction_L2_error_0003.dat' using 1:6 title "total" with lines lt 6

set title 'Spatial Reconstruction L2 Error Using 4 Modes'
plot '../results/spatial_reconstruction_L2_error_0004.dat' using 1:2 title "u" with lines lt 2,	\
	'../results/spatial_reconstruction_L2_error_0004.dat' using 1:3 title "v" with lines lt 3, \
	'../results/spatial_reconstruction_L2_error_0004.dat' using 1:4 title "w" with lines lt 4, \
	'../results/spatial_reconstruction_L2_error_0004.dat' using 1:5 title "p" with lines lt 5, \
	'../results/spatial_reconstruction_L2_error_0004.dat' using 1:6 title "total" with lines lt 6

set title 'Spatial Reconstruction L2 Error Using 5 Modes'
plot '../results/spatial_reconstruction_L2_error_0005.dat' using 1:2 title "u" with lines lt 2,	\
	'../results/spatial_reconstruction_L2_error_0005.dat' using 1:3 title "v" with lines lt 3, \
	'../results/spatial_reconstruction_L2_error_0005.dat' using 1:4 title "w" with lines lt 4, \
	'../results/spatial_reconstruction_L2_error_0005.dat' using 1:5 title "p" with lines lt 5, \
	'../results/spatial_reconstruction_L2_error_0005.dat' using 1:6 title "total" with lines lt 6

set title 'Spatial Reconstruction L2 Error Using 6 Modes'
plot '../results/spatial_reconstruction_L2_error_0006.dat' using 1:2 title "u" with lines lt 2,	\
	'../results/spatial_reconstruction_L2_error_0006.dat' using 1:3 title "v" with lines lt 3, \
	'../results/spatial_reconstruction_L2_error_0006.dat' using 1:4 title "w" with lines lt 4, \
	'../results/spatial_reconstruction_L2_error_0006.dat' using 1:5 title "p" with lines lt 5, \
	'../results/spatial_reconstruction_L2_error_0006.dat' using 1:6 title "total" with lines lt 6

set title 'Spatial Reconstruction L2 Error Using 7 Modes'
plot '../results/spatial_reconstruction_L2_error_0007.dat' using 1:2 title "u" with lines lt 2,	\
	'../results/spatial_reconstruction_L2_error_0007.dat' using 1:3 title "v" with lines lt 3, \
	'../results/spatial_reconstruction_L2_error_0007.dat' using 1:4 title "w" with lines lt 4, \
	'../results/spatial_reconstruction_L2_error_0007.dat' using 1:5 title "p" with lines lt 5, \
	'../results/spatial_reconstruction_L2_error_0007.dat' using 1:6 title "total" with lines lt 6

set title 'Spatial Reconstruction L2 Error Using 8 Modes'
plot '../results/spatial_reconstruction_L2_error_0008.dat' using 1:2 title "u" with lines lt 2,	\
	'../results/spatial_reconstruction_L2_error_0008.dat' using 1:3 title "v" with lines lt 3, \
	'../results/spatial_reconstruction_L2_error_0008.dat' using 1:4 title "w" with lines lt 4, \
	'../results/spatial_reconstruction_L2_error_0008.dat' using 1:5 title "p" with lines lt 5, \
	'../results/spatial_reconstruction_L2_error_0008.dat' using 1:6 title "total" with lines lt 6

set title 'Spatial Reconstruction L2 Error Using 9 Modes'
plot '../results/spatial_reconstruction_L2_error_0009.dat' using 1:2 title "u" with lines lt 2,	\
	'../results/spatial_reconstruction_L2_error_0009.dat' using 1:3 title "v" with lines lt 3, \
	'../results/spatial_reconstruction_L2_error_0009.dat' using 1:4 title "w" with lines lt 4, \
	'../results/spatial_reconstruction_L2_error_0009.dat' using 1:5 title "p" with lines lt 5, \
	'../results/spatial_reconstruction_L2_error_0009.dat' using 1:6 title "total" with lines lt 6

set title 'Spatial Reconstruction L2 Error Using 10 Modes'
plot '../results/spatial_reconstruction_L2_error_0010.dat' using 1:2 title "u" with lines lt 2,	\
	'../results/spatial_reconstruction_L2_error_0010.dat' using 1:3 title "v" with lines lt 3, \
	'../results/spatial_reconstruction_L2_error_0010.dat' using 1:4 title "w" with lines lt 4, \
	'../results/spatial_reconstruction_L2_error_0010.dat' using 1:5 title "p" with lines lt 5, \
	'../results/spatial_reconstruction_L2_error_0010.dat' using 1:6 title "total" with lines lt 6


# -------------------------------- EOF -------------------------------- 
