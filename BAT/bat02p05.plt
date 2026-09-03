
#set terminal wxt
#set style circle radius 10.0
#set style fill solid
set logscale x
set lmargin 10
set key off

set xrange [0.001:0.9]
set xtics (0.001,0.003,0.01,0.03,0.1,0.3,0.9)

set multiplot layout 2,1
plot 'bt02vx100x50.dat' using 1:2 with lines lc 'black',\
     'bt02vx200x50.dat' using 1:2 with lines lc 'red'
set yrange [-0.0002:0.0002]
set ytics -0.0002,0.0001,0.0002
set arrow from 0.001,0.0001 to 0.9,0.0001 nohead front lc 'grey'
set arrow from 0.001,-0.0001 to 0.9,-0.0001 nohead front lc 'grey'
plot 'bt02vx100x50.dat' using 1:5 with lines lc 'black',\
     'bt02vx200x50.dat' using 1:5 with lines lc 'red'
unset multiplot

#pause -1 "Hit any key to continue "

