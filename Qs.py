#!/usr/bin/python
from numpy import *
from scipy.interpolate import UnivariateSpline



from pylab import *
import sys

ifile=sys.argv[1]



data=loadtxt(ifile)
r=data[:,0]
D=data[:,1]
shiftD=D-(1.0-exp(-0.5))
fD=UnivariateSpline(r, shiftD, s=0.0001)



Qs=1.0/(fD.roots()[0]/sqrt(2.0))



print Qs


#f=InterpolatedUnivariateSpline(x, y)
#print T, (f.roots()[1] ),  (f.roots()[0] )  
