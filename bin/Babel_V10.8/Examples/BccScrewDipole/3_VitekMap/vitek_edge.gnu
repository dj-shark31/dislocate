# Multiplying factor
factor=20.

# Threshold in A for displacement
threshold = 0.

xScaling(x,y) = sqrt(x**2+y**2)>threshold ? x*factor : 1/0
yScaling(x,y) = sqrt(x**2+y**2)>threshold ? y*factor : 1/0

set pointsize 1
set border 0
unset ytics
unset xtics
lx =    41.8863    *1.01
ly =    16.1220    *1.01
x0 =    0.00000    -0.005*lx
y0 =    0.00000    -0.005*ly
set xrange [x0:x0+lx]
set yrange [y0:y0+ly]
set size ratio -1
 
unset key
plot 'vitek.points' u 1:2 every :3::0 w p pt 31 lc rgb 'black',\
     'vitek.points' u 1:2 every :3::1 w p pt 65 lc rgb 'black',\
     'vitek.points' u 1:2 every :3::2 w p pt 31 lc rgb 'grey',\
     'vitek_edge.res' u (0.5*($1+$4)-0.5*factor*$7):(0.5*($2+$5)-0.5*factor*$8):(xScaling($7,$8)):(yScaling($7,$8)) w vectors filled head lt 2
 
set output 'vitek_edge.eps'
set term post eps enh color solid 'Times-Roman' 20
replot
set output
set term pop
 
