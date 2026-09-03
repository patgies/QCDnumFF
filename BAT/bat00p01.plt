
#set terminal wxt
#set style circle radius 10.0
#set style fill solid
set logscale x
#set rmargin 30
set key off

set multiplot layout 2,3
plot 'batune00.dat' using 1:3 with lines lc 'black'
plot 'batune00.dat' using 1:5 with lines lc 'black'
plot 'batune00.dat' using 1:7 with lines lc 'black'
plot 'batune00.dat' using 1:4 with lines lc 'black'
plot 'batune00.dat' using 1:6 with lines lc 'black'
plot 'batune00.dat' using 1:8 with lines lc 'black'
unset multiplot

#pause -1 "Hit any key to continue "

