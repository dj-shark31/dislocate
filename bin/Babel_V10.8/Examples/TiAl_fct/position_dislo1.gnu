#!/usr/bin/gnuplot

# ratio c/a of L10 structure: 4.076/3.999 = 1.019
coa = 1.019

# Periodicity vectors in the rotated reference frame
U1x = 1./2.*sqrt(2.+4.*coa**2)
U1y = 0.
U1z = 0.

U2x = (coa**2-1.)*sqrt(2./(1.+2.*coa**2))
U2y = 3.*coa/sqrt(1.+2.*coa**2)
U2z = 0.

U3x = 0.
U3y = 0.
U3z = sqrt(2.)/2.

# Reduced coordinates of both dislocations
s1_D1 = 2.
s2_D1 = 3./2.

s1_D2 = 6.
s2_D2 = 9./2.

# Cartesian coordinates of both dislocations
x_D1 = s1_D1*U1x + s2_D1*U2x
y_D1 = s1_D1*U1y + s2_D1*U2y

x_D2 = s1_D2*U1x + s2_D2*U2x
y_D2 = s1_D2*U1y + s2_D2*U2y

show var

