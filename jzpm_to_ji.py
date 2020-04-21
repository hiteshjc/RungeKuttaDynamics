import numpy as N
import sys

which=sys.argv[1]
sqrt2=N.sqrt(2.0)

if (which=="savary"):
        Jppmm=0.042
        Jzz=-0.025
        Jpm=0.065
        Jzpm=-0.0088
        gxy=5.97
        gz=2.45

if (which=="singh"):
        # Singh
        Jppmm=0.05
        Jpm=0.05
        Jzpm=-0.14
        Jzz=0.166
        gxy=4.32
        gz=1.80


if (which=="ross"):
        # Ross
        Jppmm=0.05
        Jpm=0.05
        Jzpm=-0.14
        Jzz=0.17
        gxy=4.32
        gz=1.80

if (which=="robert"):
        # Robert
        Jppmm=0.04
        Jpm=0.085
        Jzpm=-0.15
        Jzz=0.07
        gxy=4.09
        gz=2.06

if (which=="radu"):
        # Radu Coldea
        Jppmm=0.048
        Jpm=0.074
        Jzpm=-0.159
        Jzz=0.026
        gxy=4.17
        gz=2.14

if (which=="nev"):
        # Andriy Ce
        Jppmm=-0.011
        Jpm=0.0039
        Jzpm=0.0832
        Jzz=-0.051
        gxy=4.17
        gz=2.14

J1=(1.0/3.0)*( (2.0*Jppmm)  + (4.0*Jpm) + (2.0*sqrt2*Jzpm) - (Jzz))
J2=(1.0/3.0)*( (4.0*Jppmm)  - (4.0*Jpm) + (4.0*sqrt2*Jzpm) + (Jzz))
J3=(1.0/3.0)*( (-4.0*Jppmm) - (2.0*Jpm) + (2.0*sqrt2*Jzpm) - (Jzz))
J4=(1.0/3.0)*( (2.0*Jppmm)  - (2.0*Jpm) - (1.0*sqrt2*Jzpm) - (Jzz))

print "J1 "+which+" =",J1
print "J2 "+which+" =",J2
print "J3 "+which+" =",J3
print "J4 "+which+" =",J4


