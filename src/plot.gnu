#!/usr/bin/gnuplot
set term pngcairo size 640,480 enhanced  font "Arial,12" 

reset 
load 'src/bentcoolwarm.pal'
xsize = 800
ysize = 200

set xlabel r"x"
set xrange[0:xsize]
set yrange[0:ysize]
set ylabel r"y"
set cbrange[-.04:.04]
set size ratio -1
set pm3d interpolate 2,2

!set palette viridis

i = 0
while(1){
    file = sprintf("data/w_%i.dat",i)
    set title r"w_3"
    plot file binary array=800x200 format="%float" with image notitle
    refresh
    i = i+1
    
}



