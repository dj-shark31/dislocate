# Lattice parameter
a0 =    2.94004631000       # A

# Burgers vector
b = a0

# Factor for hcp structure
factor=2./b

# Threshold in Å for displacement
threshold = 0.1*b

# Function to determine height of the atom along z 
# 	(result: 0 or 1)
cote(z) =  floor( 2.*z/a0 + 1.25 )

set border 0
unset ytics
unset xtics
set size ratio -1

set xrange [:]
set yrange [:]

show var
set size ratio -1


unset key
set para
set view map

# Color map definition for nye tensor
set cbrange [-0.25:0.25]	# (normalized density)
set palette defined ( 0 "blue", 1 "white", 2 "#DC143C" )
unset colorbox


set pointsize 2

# Style of the arrows used to plot differential displacement map
set style arrow 1 filled head lt 1 lc rgb 'black' lw 1

# Symbols used for the atoms depending on their environment and their height
set style line 11 pt 31 ps 1.5 lc rgb '#00008B' lw 1	# hcp, z=0
set style line 12 pt 65 ps 1.5 lc rgb '#00008B' lw 1	# hcp, z=1
set style line 21 pt  9 ps 2.0 lc rgb '#00008B' lw 1	# (10-11) pyramidal fault, z=0
set style line 22 pt 66 ps 2.0 lc rgb '#00008B' lw 1	# (10-11) pyramidal fault, z=1
set style line 31 pt 11 ps 2.0 lc rgb '#00008B' lw 1	# (-1011) pyramidal fault, z=0
set style line 32 pt 67 ps 2.0 lc rgb '#00008B' lw 1	# (-1011) pyramidal fault, z=1
set style line 41 pt 47 ps 1.5 lc rgb '#00008B' lw 1	# prismatic fault, z=0
set style line 42 pt 64 ps 1.5 lc rgb '#00008B' lw 1	# prismatic fault, z=1


splot 'ti_dislo_nye.fitFourier' u 1:2:(-a0*$3) w image,\
      'ti_dislo_nye.res' u ( (($4==1)||($4==2))&&(cote($3)==0)?$1:1/0):2:(0.) w p ls 11,\
      'ti_dislo_nye.res' u ( (($4==1)||($4==2))&&(cote($3)==1)?$1:1/0):2:(0.) w p ls 12,\
      'ti_dislo_nye.res' u ( (($4==3)||($4==4))&&(cote($3)==0)?$1:1/0):2:(0.) w p ls 21,\
      'ti_dislo_nye.res' u ( (($4==3)||($4==4))&&(cote($3)==1)?$1:1/0):2:(0.) w p ls 22,\
      'ti_dislo_nye.res' u ( (($4==5)||($4==6))&&(cote($3)==0)?$1:1/0):2:(0.) w p ls 31,\
      'ti_dislo_nye.res' u ( (($4==5)||($4==6))&&(cote($3)==1)?$1:1/0):2:(0.) w p ls 32,\
      'ti_dislo_nye.res' u ( (($4==7)||($4==8))&&(cote($3)==0)?$1:1/0):2:(0.) w p ls 41,\
      'ti_dislo_nye.res' u ( (($4==7)||($4==8))&&(cote($3)==1)?$1:1/0):2:(0.) w p ls 42,\
      'vitek.res' u (0.5*($1+$4)-0.5*factor*$7*($4-$1)):(0.5*($2+$5)-0.5*factor*$7*($5-$2)):(0.):(factor*$7*($4-$1)):((abs($7)>threshold)?(factor*$7*($5-$2)):1/0):(0.) w vectors arrowstyle 1


set output 'nye.eps'
set term post eps color solid 'Times-Roman' 28
replot
set output
set term pop
 
