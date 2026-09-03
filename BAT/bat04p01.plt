
#set terminal wxt
#set style circle radius 10.0
#set style fill solid
set logscale x
set logscale y
set lmargin 10
set key off

set multiplot layout 2,1
plot 'bt04xvsx1.dat' using 1:3 with lines lw 2 lc 'black',\
     'bt04xvsx2.dat' using 1:3 with lines lw 2 lc 'blue',\
     'bt04xvsx2.dat' using 1:2 with lines lw 1 lc 'grey',\
     'bt04xvsx3.dat' using 1:3 with lines lw 2 lc 'green',\
     'bt04xvsx3.dat' using 1:2 with lines lw 1 lc 'grey'
plot 'bt04xvsq1.dat' using 1:3 with lines lw 2 lc 'black',\
     'bt04xvsq1.dat' using 1:2 with lines lw 1 lc 'grey',\
     'bt04xvsq2.dat' using 1:3 with lines lw 2 lc 'blue',\
     'bt04xvsq2.dat' using 1:2 with lines lw 1 lc 'grey',\
     'bt04xvsq3.dat' using 1:3 with lines lw 2 lc 'green'
unset multiplot

#pause -1 "Hit any key to continue "

