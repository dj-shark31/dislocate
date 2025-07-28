# Lattice parameters
a0  = 2.936	# Å
coa = 1.583
c0  = coa*a0 

# Burgers vector
b = a0

# Factor for hcp structure
factor=2./b

# Threshold in A for displacement
threshold = 0.1*b

set border 0
unset ytics
unset xtics
set size ratio -1

# Positions d'origine de la dislo
x0Dislo1 = a0*sqrt(3.)*(1.+1./6.)  ; y0Dislo1 = (5.+1./4.)*c0
x0Dislo2 = a0*sqrt(3.)*(4.+1./6.)  ; y0Dislo2 = (9.+1./4.)*c0

# Positions actuelles de la dislo
xDislo1 = x0Dislo1 + a0*sqrt(3.)/4.
yDislo1 = y0Dislo1 + c0/2.

# Bords de la boîte
xmin = x0Dislo1-a0*sqrt(3.)
xmax = x0Dislo1+1.5*a0*sqrt(3.)
ymin = y0Dislo1-2.5*c0
ymax = y0Dislo1+2.5*c0
set xrange [xmin:xmax]
set yrange [ymin:ymax]

unset key
set para
set view map

# Palette pour la densité de vecteur de Burgers
#set pm3d explicit
#set cbrange [-0.1:0.1]	# non normalisé
set cbrange [-0.25:0.25]	# normalisé
set palette defined ( 0 "blue", 1 "white", 2 "#DC143C" )
unset colorbox
#set cbrange [0.:0.1]
#set palette defined ( 0 "white", 2 "#DC143C" )

# Style des flèches pour la carte de Vitek
set style arrow 1 filled head lt 1 lc rgb 'black' lw 2

# Symboles pour les atomes
set style line 11 pt 31 ps 2 lc rgb '#00008B' lw 2	# hcp, z=0
set style line 12 pt 65 ps 2 lc rgb '#00008B' lw 2	# hcp, z=1
set style line 21 pt 63 ps 2 lc rgb '#00008B' lw 2	# twin 1, z=0
set style line 22 pt 68 ps 3 lc rgb '#00008B' lw 2	# twin 1, z=1
set style line 31 pt  2 ps 2 lc rgb '#00008B' lw 5	# twin 2, z=0 (absent normalement)
set style line 32 pt  3 ps 2 lc rgb '#00008B' lw 5	# twin 2, z=1
set style line 41 pt 47 ps 2 lc rgb '#00008B' lw 2	# faute prism., z=0
set style line 42 pt 64 ps 2 lc rgb '#00008B' lw 2	# faute prism., z=1

# Symboles pour le centre de la dislocation
set style line 1 pt 1 lc rgb '#DC143C' ps 3 lw 3

# Cote des atomes suivant z (0 ou 1)
cote(z) =  floor( 2.*z/a0 + 0.25 )

splot "ti_dislo_screw.fitFourier" u 1:2:(-a0*$3) w pm3d,\
      xDislo1,yDislo1,0 t '' w p ls 1,\
      'ti_dislo_screw.res' u 1:2:( ( ($4==1)||($4==2) )&&(cote($3)==0)?0:1/0 ) t '' w p ls 11,\
      'ti_dislo_screw.res' u 1:2:( ( ($4==1)||($4==2) )&&(cote($3)==1)?0:1/0 ) t '' w p ls 12,\
      'ti_dislo_screw.res' u 1:2:( ( ($4==3)||($4==4) )&&(cote($3)==0)?0:1/0 ) t '' w p ls 21,\
      'ti_dislo_screw.res' u 1:2:( ( ($4==3)||($4==4) )&&(cote($3)==1)?0:1/0 ) t '' w p ls 22,\
      'ti_dislo_screw.res' u 1:2:( ( ($4==5)||($4==6) )&&(cote($3)==0)?0:1/0 ) t '' w p ls 31,\
      'ti_dislo_screw.res' u 1:2:( ( ($4==5)||($4==6) )&&(cote($3)==1)?0:1/0 ) t '' w p ls 32,\
      'ti_dislo_screw.res' u 1:2:( ( ($4==7)||($4==8) )&&(cote($3)==0)?0:1/0 ) t '' w p ls 41,\
      'ti_dislo_screw.res' u 1:2:( ( ($4==7)||($4==8) )&&(cote($3)==1)?0:1/0 ) t '' w p ls 42,\
      'vitek.res' u (0.5*($1+$4)-0.5*factor*$7*($4-$1)):(0.5*($2+$5)-0.5*factor*$7*($5-$2)):(2.):(factor*$7*($4-$1)):((abs($7)>threshold)?(factor*$7*($5-$2)):1/0):(0.) w vectors arrowstyle 1
 
set output 'vitek_screw.eps'
set term post eps enh color solid 'Times-Roman' 20
replot
set output
set term pop
 
show var
