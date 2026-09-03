
#set terminal wxt
#set style circle radius 10.0
#set style fill solid
set logscale x
set lmargin 10
set key off

#set yrange [-0.0002:0.0002]
#set ytics -0.0002,0.0001,0.0002

set multiplot layout 2,1
set yrange [0.0:0.8]
set ytics out
set ytics 0,0.2,0.8
set xrange [0.1:0.9]
unset xtics
set tmargin at screen 0.95
set bmargin at screen 0.68
plot 'bt02ghx100x50.dat' using 1:2 with lines lc 'black',\
     'bt02ghx200x50.dat' using 1:2 with lines lc 'red'
set tmargin at screen 0.68
set bmargin at screen 0.55
set yrange [-1.5e-4:5e-5]
set ytics -1.5e-4,1.5e-4,3e-5
set xtics out
set xtics (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
plot 'bt02ghx100x50.dat' using 1:4 with lines lc 'black',\
     'bt02ghx200x50.dat' using 1:4 with lines lc 'red'
     
set yrange [0.0:1.5]
set ytics out
set ytics 0,0.5,1.5
unset xtics
set tmargin at screen 0.48
set bmargin at screen 0.21
plot 'bt02shx100x50.dat' using 1:2 with lines lc 'black',\
     'bt02shx200x50.dat' using 1:2 with lines lc 'red'
set tmargin at screen 0.21
set bmargin at screen 0.08
set xtics out
set xtics (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
set yrange [-1.5e-4:5e-5]
set ytics -1.5e-4,1.5e-4,3e-5
plot 'bt02shx100x50.dat' using 1:4 with lines lc 'black',\
     'bt02shx200x50.dat' using 1:4 with lines lc 'red'
unset multiplot

#pause -1 "Hit any key to continue "

