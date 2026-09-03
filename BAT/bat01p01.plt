
#set terminal wxt
#set style circle radius 10.0
#set style fill solid
set logscale x
#set rmargin 30
set key off

set multiplot layout 2,2
plot 'bt01vsx1.dat' using 1:2 with lines lc 'red' ,\
     'bt01vsx2.dat' using 1:2 with lines lc 'black'
plot 'bt01vsx1.dat' using 1:3 with lines lc 'red' ,\
     'bt01vsx2.dat' using 1:3 with lines lc 'black'
plot 'bt01vsq1.dat' using 1:2 with lines lc 'red' ,\
     'bt01vsq2.dat' using 1:2 with lines lc 'blue',\
     'bt01vsq3.dat' using 1:2 with lines lc 'black'
plot 'bt01vsq1.dat' using 1:3 with lines lc 'red' ,\
     'bt01vsq2.dat' using 1:3 with lines lc 'blue',\
     'bt01vsq3.dat' using 1:3 with lines lc 'black'
unset multiplot

#pause -1 "Hit any key to continue "

