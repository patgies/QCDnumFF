
#set terminal wxt
#set style circle radius 10.0
#set style fill solid
set logscale x
set lmargin 10
set rmargin 10
set xtics (0.001,0.003,0.01,0.03,0.1,0.3,1)
set key off

set multiplot layout 3,1
set yrange [0:3]
set ytics 0,1,3
set tmargin at screen 0.95
set bmargin at screen 0.70
plot 'batune00.dat' using 1:13 with lines lc 'black',\
     'batune00.dat' using 1:15 with lines lc 'red'
set tmargin at screen 0.63
set bmargin at screen 0.38
plot 'batune00.dat' using 1:13 with lines lc 'black',\
     'batune00.dat' using 1:16 with lines lc 'red'
set yrange [0:1]
set ytics 0,0.2,1
set tmargin at screen 0.31
set bmargin at screen 0.06
plot 'batune00.dat' using 1:14 with lines lc 'black',\
     'batune00.dat' using 1:17 with lines lc 'red'
unset multiplot

#pause -1 "Hit any key to continue "

