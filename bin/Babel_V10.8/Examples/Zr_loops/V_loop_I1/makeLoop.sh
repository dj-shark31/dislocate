#!/bin/bash

# Exécutables
BABEL="$HOME/Babel/bin/babel"

# Paramètre de maille et rapport c/a
alat=3.234
coa=1.598021026592455164

# Vecteur de Burgers
bx="0."
by=".57735026918962576451"
bz=`echo "-$coa/2." | bc -l`

# Taille de la boucle
n=5

# Fichier avec la structure hcp parfaite
hcpFile="zr_duplicated.xyz"

# Nombre d'atomes
imm=`head -1 $hcpFile`

# Calcul des 6 extrémités de la boucle
xLoop1=`echo "0.5*$n+0.25" | bc -l`
yLoop1=`echo "sqrt(3.)/2.*($n+1./6.)" | bc -l`

xLoop2=`echo "$n+0.25" | bc -l`
yLoop2=`echo "sqrt(3.)/2.*1./6." | bc -l`

xLoop3=`echo "0.5*$n" | bc -l`
yLoop3=`echo "-sqrt(3.)/2.*($n+1./3.)" | bc -l`

xLoop4=`echo "-0.5*$n" | bc -l`
yLoop4=`echo "-sqrt(3.)/2.*($n+1./3.)" | bc -l`

xLoop5=`echo "-$n-0.25" | bc -l`
yLoop5=`echo "sqrt(3.)/2.*1./6." | bc -l`

xLoop6=`echo "-0.5*$n-0.25" | bc -l`
yLoop6=`echo "sqrt(3.)/2.*($n+1./6.)" | bc -l`

zLoop=0.

# Création de la boucle
echo "
&input
  alat=$alat
  imm=$imm
  inpXyz=.true.
  inpFile="$hcpFile"
  out_alat=$alat
  outCfg=.true.
  outFile="zr_loop.cfg"

  remove_cut=.true.

  initial=.false.

  xImages=.true.
  yImages=.true.
  zImages=.true.
  nxImages=1
  nyImages=1
  nzImages=1

  !max_Euler=10

  verbosity=5
//

&loops

  ! Noise for loop points
  xLoop_noise=0.d-6

  ! Number of loops
  nLoop=1

  ! Burgers vector of loop 1
  bLoop(1,1)=$bx
  bLoop(2,1)=$by
  bLoop(3,1)=$bz

  ! Number of points in loop 1
  iLoop(1)=6

  ! Coordinates of points for loop 1
  xLoop(1,1,1)=$xLoop1
  xLoop(2,1,1)=$yLoop1
  xLoop(3,1,1)=$zLoop

  xLoop(1,2,1)=$xLoop2
  xLoop(2,2,1)=$yLoop2
  xLoop(3,2,1)=$zLoop

  xLoop(1,3,1)=$xLoop3
  xLoop(2,3,1)=$yLoop3
  xLoop(3,3,1)=$zLoop

  xLoop(1,4,1)=$xLoop4
  xLoop(2,4,1)=$yLoop4
  xLoop(3,4,1)=$zLoop

  xLoop(1,5,1)=$xLoop5
  xLoop(2,5,1)=$yLoop5
  xLoop(3,5,1)=$zLoop

  xLoop(1,6,1)=$xLoop6
  xLoop(2,6,1)=$yLoop6
  xLoop(3,6,1)=$zLoop

//
" > input_loop.babel
$BABEL input_loop.babel | tee output_loop.babel
