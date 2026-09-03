
#set terminal wxt
#set style circle radius 10.0
#set style fill solid
set logscale x
#set rmargin 30
set key off
set yrange [-0.001:0.001]
set ytics out
set ytics -0.001,0.0005,0.001
set xtics out
set xrange [0.1:0.9]
set xtics (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

set multiplot layout 2,1
set arrow from 0.1,0.0001 to 0.9,0.0001 nohead front lc 'grey'
set arrow from 0.1,-0.0001 to 0.9,-0.0001 nohead front lc 'grey'
plot 'bt02ghx100x50.dat' using 1:5 with lines lc 'black',\
     'bt02ghx200x50.dat' using 1:5 with lines lc 'red'
plot 'bt02shx100x50.dat' using 1:5 with lines lc 'black',\
     'bt02shx200x50.dat' using 1:5 with lines lc 'red'
unset multiplot

#pause -1 "Hit any key to continue "

