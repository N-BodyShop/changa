#gnuplot script


set terminal png #fontfile "/usr/share/texmf-tetex/fonts/type1/urw/helvetic/uhvr8a.pfb" 14 
set size 1.0,0.5
set output 'scatter.png'

# set title "Execution Time Scaling"
set title "Distance vs. Communication Volume"
set xlabel "Distance"
set ylabel "Nodes exchanged"
#set y2label "Data transferred/request"
#unset key
#set key left
set autoscale
#unset logscale; 
set ytics nomirror
#set logscale x 2
#set logscale y 2
set format x "%1.4f"
#set label "Optimal size = 32.51" at first 9.51, first 0.4822
#set arrow from 18.5, .482 to 32.51,0.7822 

set style line 1 lt rgb "red" lw 3 pt 2
set style line 2 lt rgb "blue" lw 3 pt 4
#set style line 3 lt rgb "orange" lw 3 pt 6
set style line 3 lt rgb "green" lw 3 pt 6

plot [] [] \
      'scatter.dat' using 1:2 notitle with points  ls 1
      #,\
      #'bucket_kernel.dat' using 2:5 title 'Data per request' axes x1y2 with linespoints  lw 3
