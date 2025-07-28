set style textbox opaque noborder

# Lattice parameter
a0 =    2.85000000000       # A

# Burgers vector
b = sqrt(3.)/2.*a0

# Factor for bcc structure
factor=3./b

# Threshold in A for displacement
threshold = 0.1*b

scaling(x) = abs(x)>threshold ? x*factor : 1/0

# Distances between atomic columns
lx = sqrt(6.)/3.*a0 ; ly = sqrt(2.)/2.*a0

# Box origin
x0 = 0.*lx ; y0 = 0.*ly

# Periodicity vectors
u1x = 7.5*lx ; u1y =-9.*ly
u2x = 7.5*lx ; u2y = 9.*ly

# Dislocation position
xD1 = 5.5*lx  ; yD1 = (4.+1./3.)*ly 
xD2 = xD1+0.5*(u1x-u2x) ; yD2 = yD1+0.5*(u1y-u2y)+ly/3.


set pointsize 1
set border 0
unset xtics
unset ytics
xMin = x0 
xMax = x0 + u1x + u2x
yMin = y0 + u1y
yMax = y0 + u2y
show var
set xrange [xMin:xMax]
set yrange [yMin:yMax]
set size ratio -1
 
unset key
set cbrange [-0.30/a0:0.30/a0]
set palette defined ( 0 "blue", 1 "white", 2 "#DC143C" )
unset colorbox

set para

# Atom symbols
set style line 11 pt 31 lc rgb 'black'
set style line 12 pt 65 lc rgb 'black'
set style line 13 pt 31 lc rgb 'grey' 

# Dislocation symbols
set style line 1 pt 1 lc rgb 'black' ps 1 lw 2
set style line 2 pt 1 lc rgb 'black' ps 1 lw 2


plot "fe_dislo_nye.fit" u 1:2:3 w image,\
     "fe_dislo_nye.res" u 1:2 w p ls 12 ps 2,\
     xD1,yD1         w p ls 1,\
     xD2,yD2         w p ls 2
 
set output 'nye.eps'
set term post eps enh dl 1 color 'Times-Roman' 24
replot
set output
set term pop
 
