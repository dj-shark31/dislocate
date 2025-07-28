#!/usr/bin/gnuplot

# ratio c/a of L10 structure: 4.076/3.999 = 1.019
coa = 1.019

# Components of rotation matrix
rot11 = 1./sqrt(2.+4.*coa**2) * -1.
rot12 = 1./sqrt(2.+4.*coa**2) *  1.
rot13 = 1./sqrt(2.+4.*coa**2) *  2.*coa

rot21 = 1./sqrt(1.+2.*coa**2) *  coa
rot22 = 1./sqrt(1.+2.*coa**2) * -coa
rot23 = 1./sqrt(1.+2.*coa**2) *  1.

rot31 = 1./sqrt(2.)
rot32 = 1./sqrt(2.)
rot33 = 0.

show var 
