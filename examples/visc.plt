set size square
set key right bottom
set key Left reverse
set key box
set xlabel 'Time'
set ylabel 'Uniaxial elongational viscosity'
set title 'N=100, M=100, uniaxial strain-rate = 0.001'
erate=0.001
dt=0.01


p 'press_sg.txt' u (dt*$1):($2/erate) t 'Raw data' pt 65 lc rgb "blue",\
'press_sg.txt' u (dt*$1):($3/erate) t 'Smoothed data' pt 66 lc rgb "red"
pause -1

set terminal postscript enhanced color eps ",20"
set output 'visc.eps'
replot
!convert -density 300 visc.eps visc.png
!eog visc.png

