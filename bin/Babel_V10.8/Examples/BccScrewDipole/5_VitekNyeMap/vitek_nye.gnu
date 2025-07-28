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
x0 = -6.0*lx ; y0 = 13.5*ly

# Periodicity vectors
u1x = 7.5*lx ; u1y =-9.*ly
u2x = 7.5*lx ; u2y = 9.*ly
set arrow from x0, y0 to x0+u1x, y0+u1y filled head front 
set arrow from x0, y0 to x0+u2x, y0+u2y filled head front
set arrow from x0+u1x, y0+u1y to x0+u1x+u2x, y0+u1y+u2y dt 2 lw 2 nohead front
set arrow from x0+u2x, y0+u2y to x0+u1x+u2x, y0+u1y+u2y dt 2 lw 2 nohead front
set label "U_1" at x0+u1x, y0+u1y offset 0., 1. center box front
set label "U_2" at x0+u2x, y0+u2y offset -1.,0.5 right box front

# Dislocation position
xD1 = 5.5*lx - u1x ; yD1 = (4.+1./3.)*ly - u1y
xD2 = xD1+0.5*(u1x+u2x) ; yD2 = yD1+0.5*(u1y+u2y)+ly/3.
set label "+b" at xD1, yD1 center box offset 0., -2. front
set label "-b" at xD2, yD2 left   box offset 0., -1.5 front

# Dislocation cut
#    center
xC = 0.5*(xD1+xD2) ; yC = 0.5*(yD1+yD2)
#    direction
Ax = yD1-yD2 ; Ay=xD2-xD1
set arrow from xD1, yD1 to xD2, yD2 nohead front lw 4
set arrow from xD1, yD1 to xD2, yD2 nohead front lc rgb 'white' lw 0.8
set arrow from xC, yC to xC+Ax, yC+Ay head filled front
set label "A" at xC+0.8*Ax, yC+0.8*Ay offset 0.7, 0. left box front

# Dipole vector
dy=-0.5*ly
set arrow from xD1, yD1+dy to xD2, yD2+dy filled head front lc "white" lw 2
set arrow from xD1, yD1+dy to xD2, yD2+dy filled head front 
set label "d" at 0.5*(xD1+xD2), 0.5*(yD1+yD2)+dy offset 0.,-0.5 center box front

set pointsize 1
set border 0
unset xtics
unset ytics
xMin = x0 - 6.5*ly
xMax = x0 + 17.5*ly
yMin = y0 - 12.*ly
yMax = y0 + 12.*ly
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

# Arrow style
set style arrow 1 filled head lt 1 lc rgb 'purple'

# Dislocation symbols
set style line 1 pt 1 lc rgb 'black' ps 1 lw 2
set style line 2 pt 1 lc rgb 'black' ps 1 lw 2

# Offset to plot Vitek map
dx=u1x+u2x
dy=u1y+u2y

plot "fe_dislo_nye.fit" u 1:2:3 w image,\
     'vitek.points' u ($1+dx):($2+dy) every :3::0 w p ls 11,\
     'vitek.points' u ($1+dx):($2+dy) every :3::1 w p ls 12,\
     'vitek.points' u ($1+dx):($2+dy) every :3::2 w p ls 13,\
     'vitek.res' u (0.5*($1+$4)-0.5*factor*$7*($4-$1)+dx):(0.5*($2+$5)-0.5*factor*$7*($5-$2)+dy):(scaling($7)*($4-$1)):(scaling($7)*($5-$2)) w vectors arrowstyle 1,\
     xD1,yD1         w p ls 1,\
     xD1+u1x,yD1+u1y w p ls 1,\
     xD1-u1x,yD1-u1y w p ls 1,\
     xD1+u2x,yD1+u2y w p ls 1,\
     xD1-u2x,yD1-u2y w p ls 1,\
     xD2,yD2         w p ls 2,\
     xD2-u1x,yD2-u1y w p ls 2,\
     xD2-u2x,yD2-u2y w p ls 2,\
     xD2-u1x-u2x,yD2-u1y-u2y w p ls 2
 
set output 'vitek_nye.eps'
set term post eps enh dl 1 color 'Times-Roman' 24
replot
set output
set term pop
 
