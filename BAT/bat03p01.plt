
#set terminal wxt
#set style circle radius 10.0
#set style fill solid
set logscale x
#set logscale y
set lmargin 20
set rmargin 20
set key off

set multiplot layout 3,1
set yrange [0:0.6]
set ytics 0,0.2,0.6
set xrange[0.001:0.9]
set xtics (0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.9)
plot 'bt03FLx1.dat' using 1:2 with lines lc 'blue',\
     'bt03FLx3.dat' using 1:2 with lines lc 'red'
set yrange [0:5]
set ytics 0,1,5
plot 'bt03F2x1.dat' using 1:2 with lines lc 'blue',\
     'bt03F2x3.dat' using 1:2 with lines lc 'red'
set yrange [0:0.8]
set ytics 0,0.2,0.8
plot 'bt03F3x1.dat' using 1:2 with lines lc 'blue',\
     'bt03F3x3.dat' using 1:2 with lines lc 'red'
unset multiplot

#pause -1 "Hit any key to continue "

