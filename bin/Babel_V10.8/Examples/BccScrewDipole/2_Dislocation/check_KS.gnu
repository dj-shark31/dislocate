# Elastic constants in cubic reference frame
  C11=232.27	# GPa
  C12=135.87
  C44=101.9
  Cp=(C11-C12)/2.

# Energy prefactor for screw dislo
Ks=sqrt(C44*Cp)	# GPa
Ks=(C44+2*Cp)/3.

show var Ks

