
#set terminal wxt
#set style circle radius 10.0
#set style fill solid
set logscale x
#set rmargin 30
set key off
set yrange [-0.0002:0.0002]
set ytics out
set ytics -0.0002,0.0001,0.0002
set xtics out
set xrange [0.001:0.9]
set xtics (0.001,0.003,0.01,0.03,0.1,0.3,0.9)

set multiplot layout 2,1
set arrow from 0.001,0.0001 to 0.9,0.0001 nohead front lc 'grey'
set arrow from 0.001,-0.0001 to 0.9,-0.0001 nohead front lc 'grey'
plot 'bt02gx100x50.dat' using 1:5 with lines lc 'black',\
     'bt02gx200x50.dat' using 1:5 with lines lc 'red'
plot 'bt02sx100x50.dat' using 1:5 with lines lc 'black',\
     'bt02sx200x50.dat' using 1:5 with lines lc 'red'
unset multiplot

#pause -1 "Hit any key to continue "

