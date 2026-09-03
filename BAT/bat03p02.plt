
#set terminal wxt
#set style circle radius 10.0
#set style fill solid
set logscale x
set logscale y
#set rmargin 40
set key off
set view map
set xtics out
set ytics out
set style data pm3d
set xrange[0.001:0.9]
set xtics (0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.9)

set multiplot layout 2,1
set lmargin at screen 0.10
set rmargin at screen 0.85
set bmargin at screen 0.45
set tmargin at screen 0.92
set yrange [100:30000]
splot 'bt03F3xq.dat' using 1:2:3,\
      'bt03ndxq.dat' using 1:2:3 with lines lc 2 lw 1 lt 5 ,\
      'bt03rsxq.dat' using 1:2:3 with lines lc 3 lw 3 lt 5
set bmargin at screen 0.10
set tmargin at screen 0.35
unset logscale y
set yrange [-0.0015:0.0015]
set ytics -0.0015,0.0005,0.0015
set arrow from 0.001,0.0005 to 0.9,0.0005 nohead front lc 'grey'
set arrow from 0.001,-0.0005 to 0.9,-0.0005 nohead front lc 'grey'
plot 'bt03F3x1.dat' using 1:6 with lines lc 'black',\
     'bt03F3x3.dat' using 1:6 with lines lc 'red'
unset multiplot
#pause -1 "Hit any key to continue "

