#!/bin/sh

# Needed programs
#   program displacement (to compile in ~/Babel/Source)
DISPLACEMENT="${HOME}/Babel/bin/displacement"
#   program patternInit (to compile in ~/Babel/Source)
PATTERNINIT="${HOME}/Babel/bin/patternInit"
#   program vitek (to compile in ~/Babel/Source/Tools/Vitek
VITEK="${HOME}/Babel/bin/vitek_hcp"
#   program fourierFit2D (to compile in ~/Babel/Source/Tools/FourierInterpolate)
FOURIERFIT2D="${HOME}/Babel/bin/fourierFit2D"

# Input file with the dislocations
#    right now, I use the input file which has been generated with Babel
#    but you should of course use the one obtained after relaxation with Vasp
#inpFile='out.cfg'
inpFile='POSCAR_out'

# Reference file without the dislocation
refFile='pure.cfg'

# Lattice parameters
alat="2.940046310"	# Å
coa="1.57959627589070942975"

# Maximal number of atoms
nAtoms="576"

# Dimension of the simulation cell in x and z directions
#  (needed only for Fourier interpolation of dislocation density)
lx=41.79677582
lz=40.66209412
# number of Fourier coefficients in each direction
#    it cannot be larger than the number of periodic vectors in this direction
#    it has to be decreases if the interpolation is not smooth
nx=9
nz=8

#################################################
# Step 1: the Vitek differential displacement map
#
# Calculate the displacement between the input and the reference files 
#    and also rotate the crystal because the program which draws the Vitek map
#    assumes that the dislocation line is along z axis
#
# name of the file containing displacement
disp1File='ti_dislo_displacement.xyz'
# creation of the input file for program displacement
cat > input_vitek.displacement << EOF
&input
  
  ! Rotation
  rotate=.true.
  rot(1,1)=1.  
  rot(1,2)=0.  
  rot(1,3)=0.
  rot(2,1)=0.  
  rot(2,2)=0.  
  rot(2,3)=1.
  rot(3,1)=0.  
  rot(3,2)=-1. 
  rot(3,3)=0.

  ! Input structure (with the dislocations)
  inpFile='${inpFile}'
  inpPoscar=.true.
  !inpCfg=.true.

  ! Reference structure (perfect crystal)
  refFile='${refFile}'
  refCfg=.true.

  ! Output file containing displacement
  outFile='${disp1File}'
  outXyz=.true.
  initial=.true.	
//
EOF
# launch program displacement
${DISPLACEMENT} input_vitek.displacement

# Produce all files needed to draw Vitek map
#  this will generate the files:
#     - vitek.points which contains the atoms positions gathered in layers
#     - vitek.res which contains the arrows definitions
#     - vitek.gnu: the gnuplot script used to plot the DD map
#		this script can be modified
#		if this script already exist, the vitek program won't erase it
${VITEK} ${disp1File} ${disp1File} ${alat} perio

# Draw Vitek map (not needed)
# 	this produces the image file vitek.eps
gnuplot --persist vitek.gnu


#################################################
# Step 2: the Nye tensor calculation
#
# Calculate the dislocation density (projection of the Nye tensor along the direction 
#     of the screw component)
#     The same rotation is needed
#     We also need to duplicate the structure in the dislocation direction to be able
#     to find all atom neighbours
#
# Creation of the file containing the definition of perfect atomic environment for hcp lattice
#    another rotation is needed becasue the program assumes a different orientation (see online doc)
cat > input.patternInit << EOF
&input

  crystal='hcp'
  alat=${alat}
  coa=${coa}

  ! Rotation
  rotate=.true.
  rot(1,1)=0.  
  rot(1,2)=0.  
  rot(1,3)=1.
  rot(2,1)=0.  
  rot(2,2)=1.  
  rot(2,3)=0.
  rot(3,1)=1.  
  rot(3,2)=0. 
  rot(3,3)=0.

  ! Patterns for pyramidal and prismatic faults
  addHcpPy1Fault=.true.
  addHcpPrismFault=.true.

  ! File where to store the patterns
  patternFile="hcp.pattern"

//
EOF
# lauch program patternInit
${PATTERNINIT} input.patternInit


# Number of atoms after duplication
imm=`echo "$nAtoms * 2" | bc -l`
echo $imm
#
# name of the file containing displacement
disp2File='ti_dislo_nye.res'
#
# creation of the input file for program displacement
cat > input_nye.displacement << EOF
&input
 
  ! Duplication of the structure
  !   (this will be performed before rotation)
  duplicate=.true.
  lat(1)=1
  lat(2)=2
  lat(3)=1

  ! Number of atoms after duplication
  imm=${imm}
 
  ! Rotation
  rotate=.true.
  rot(1,1)=1.  
  rot(1,2)=0.  
  rot(1,3)=0.
  rot(2,1)=0.  
  rot(2,2)=0.  
  rot(2,3)=1.
  rot(3,1)=0.  
  rot(3,2)=-1. 
  rot(3,3)=0.

  ! radius of the neighbourhood sphere
  rNeigh=1.15	! (normalized by alat)

  ! Scaling factor for the distances (lattice parameter)
  alat=${alat}
  
  ! Threshold angle for pairing vectors between input and reference structures
  patternAngleThreshold=10.d0	! degrees

  patternSelectionMethod=4

  ! Pattern is read from file
  patternFile='hcp.pattern'


  ! Input structure (with the dislocations)
  inpFile='${inpFile}'
  inpPoscar=.true.
  !inpCfg=.true.

  ! Reference structure (perfect crystal)
  refFile='${refFile}'
  refCfg=.true.

  ! Output file containing displacement
  outFile='${disp2File}'
  outOnlyAtoms=.true.
  initial=.true.	

  ! Does not output atomic displacement (default: .true.)
  out_displacement=.false.

  ! Output dislocation density corresponding to the directions defined by lNye(:) and bNye(:)
  out_BurgersDensity=.true.
  ! line direction for the dislocation density (will be normalized)
  lNye(1) = 0.d0
  lNye(2) = 0.d0
  lNye(3) = 1.d0
  ! Burgers vector direction for the dislocation density (will be normalized)
  bNye(1) = 0.d0
  bNye(2) = 0.d0
  bNye(3) = 1.d0

  ! Output closest pattern for each atom
  out_pattern=.true.

//
EOF
# launch program displacement to calculate dislocation density on each atom position
${DISPLACEMENT} input_nye.displacement


#################################################
# Step 3: interpolation of the Nye tensor with Fourier series
#   you can use a different interpolation if you want, like cubic splines (with gnuplot for instance)
#   or you can use a smearing function, like a Gaussian, to broaden the information 
#   which you only know at atomic positions
#   The purpose is to obtain a continuous 2D function, so as to make a nice drawing then!
#   The advantage of Fourier series is that you take full account of periodicity 
#   so you don't have any trouble at the borders of your simulation cell
#   But do not use too many harmonics in the Fourier series, otherwise the interpolation
#   will not be smooth anymore.
#

# create input file for 2D Fourier interpolation
cat > input_nye.fourierFit2D << EOF
${nx} ${nz}
${lx} 0.
0. ${lz}
${disp2File}
1
5
0.
1.
ti_dislo_nye.coefFourier
0. 0.
${lx} 0.
0. ${lz}
200 200
ti_dislo_nye.fitFourier
EOF
# launch program fourierFit2D to interpolate dislocation density
${FOURIERFIT2D} < input_nye.fourierFit2D

# Draw Nye tensor map with Vitek map
# 	this produces the image file nye.eps
#       (you can remove the blank using fixbb script: http://www.gnuplot.info/scripts/files/fixbb)
gnuplot --persist nye.gnu
fixbb nye.eps
