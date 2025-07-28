# Lattice parameter
a0 =    2.85000000000       # A

# Burgers vector
b = sqrt(3.)/2.*a0

# Factor for bcc structure
factor=3./b

# Threshold in A for displacement
threshold = 0.1*b

scaling(x) = abs(x)>threshold ? x*factor : 1/0

set pointsize 1
set border 0
unset ytics
unset xtics
lx =    129.149    *1.01
ly =    106.808    *1.01
x0 =   -52.3578    -0.005*lx
y0 =    0.00000    -0.005*ly
set xrange [x0:x0+lx]
set yrange [y0:y0+ly]
set size ratio -1
 
unset key
plot 'vitek.points' u 1:2 every :3::0 w p pt 31 lc rgb 'black',\
     'vitek.points' u 1:2 every :3::1 w p pt 65 lc rgb 'black',\
     'vitek.points' u 1:2 every :3::2 w p pt 31 lc rgb 'grey',\
     'vitek.res' u (0.5*($1+$4)-0.5*factor*$7*($4-$1)):(0.5*($2+$5)-0.5*factor*$7*($5-$2)):(scaling($7)*($4-$1)):(scaling($7)*($5-$2)) w vectors filled head lt 1
 
set output 'vitek.eps'
set term post eps enh color solid 'Times-Roman' 20
replot
set output
set term pop
 
