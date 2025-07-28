#!/bin/bash

BABEL="${HOME}/Babel/bin/babel"
inpFile="input.babel"
outFile="output.babel"
pi="3.14159265358979"

mag="1"



echo "#"
echo "# 1: theta (°)"
echo "# 2: pre-logarithmic factor (eV.A^-1)"
echo "# 3: core traction constant (eV.A^-1)"
echo "# 4: core dilatation contribution (eV.A)"

# Loop on angles #############
theta_min=`expr -1 \* $mag`
theta_max=`expr 361 \* $mag`
theta_mag="$theta_min"
while [ "$theta_mag" -le "$theta_max" ]; do

theta=`bc << EOF
scale=1
$theta_mag/$mag
EOF
`

cos=`bc -l << EOF
  scale=12
  c($theta*$pi/180.)
EOF
`
sin=`bc -l << EOF
  scale=12
  s($theta*$pi/180.)
EOF
`
moins_sin=`bc -l << EOF
  scale=12
  -s($theta*$pi/180.)
EOF
`

cat > $inpFile << EOF
	&input

	  anisotropic_elasticity=.true.

          ! Elastic constants (in GPa)
	  CVoigt(1,1)=266.00
          CVoigt(1,2)=140.00
          CVoigt(1,3)=134.00
          CVoigt(1,5)=-8.4853
          CVoigt(2,1)=140.00
          CVoigt(2,2)=266.00
          CVoigt(2,3)=134.00
          CVoigt(2,5)=8.4853
          CVoigt(3,1)=134.00
          CVoigt(3,2)=134.00
          CVoigt(3,3)=272.00
          CVoigt(4,4)=57.000
          CVoigt(4,6)=8.4853
          CVoigt(5,1)=-8.4853
          CVoigt(5,2)=8.4853
          CVoigt(5,5)=57.000
          CVoigt(6,4)=8.4853
          CVoigt(6,6)=63.000

          ! ! Noise to give to elastic constants to withdraw any isotropy
          ! CVoigt_noise=1.d-4	! no unit 

          ! Lattice parameter (in Å)
	  alat=2.8803

//

&dislo
	  ! Number of dislocations
	  nd=1

          ! Burgers vector (normalized by alat if defined)
	  bDislo(1,1)=0.
	  bDislo(2,1)=0.
	  bDislo(3,1)=0.8660254037844386468

          ! Line vector
	  lDislo(1,1)=${sin}
	  lDislo(2,1)=0.
	  lDislo(3,1)=${cos}

          ! Cut direction 
	  cutDislo(1,1)=${cos}
	  cutDislo(2,1)=0.
	  cutDislo(3,1)=${moins_sin}
  
          ! Dislo position
	  cDislo(1,1)=0.
	  cDislo(2,1)=0.
	  cDislo(3,1)=0.
//
EOF

$BABEL $inpFile > $outFile

Kfactor=`grep -A 2 "prelogarithmic energy factor: E = K\*ln(R/r)" $outFile | tail -1 \
	| sed -e 's/=//' -e 's/-  units = energy \/ distance (eV.A^-1)//' \
	| awk '{print $1}'`
E0=`grep -A 2 "angular dependence of the energy" $outFile | tail -1 \
	| sed -e 's/=//' -e 's/-  units = energy \/ distance (eV.A^-1)//' \
	| awk '{print $1}'`
Edipole=`grep -A 2 "energy prefactor: E = -K/r^2" $outFile | tail -1 \
	| sed -e 's/=//' -e 's/  -  units = energy * distance (eV.A)//' \
	| awk '{print $1}'`


echo -e "$theta\t$Kfactor\t$E0\t$Edipole"

theta_mag=`expr $theta_mag + 1`
done
