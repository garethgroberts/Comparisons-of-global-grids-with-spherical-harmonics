#!/bin/bash

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import re

minv = float(sys.argv[1])
maxv = float(sys.argv[2])
traw = sys.argv[3]

#t = float(re.findall(r'\d+', traw)[0])
t = float(re.findall(r"[-+]?(?:\d*\.*\d+)", traw)[0]) # must not be in scientific notation!

print(minv,maxv,t)

xyzpd = pd.read_csv('xyz.temp', sep='\t', header=None, names = ['x','y','z'])
#xyzpd = pd.read_csv('xyz.temp', sep=' ', header=None, names = ['x','y','z'])

xyzpd['latfac'] = np.power(np.power(np.cos( np.radians(xyzpd['y'])  ),2.),0.5)

#bins = np.arange(minv,maxv,t) # limited to >=minv, <=maxv
#bins = np.arange(minv-(t/2),maxv+(3*t/2),t) # centred on zero if minv, maxv symmetrical
bins = np.arange(minv+(t/2),maxv+(t/2),t) # if minv = -10, maxv = 10, these go from -9.5 to 9.5
sums = []
#bincentres = [minv-(t/2.)]
bincentres = [minv]

### below min? ### - put all into minv bin (centred)
bmin_bool = xyzpd['z'].lt(min(bins))
subset = xyzpd[bmin_bool.values]
sums.append(np.sum(subset['latfac']))

### in-between bins ###
for i in range(len(bins)-1):
	j=i+1
	v1 = bins[i]
	v2 = bins[j]
	subset = xyzpd[ (xyzpd['z'] >= v1) & (xyzpd['z'] < v2) ]
	sums.append(np.sum(subset['latfac']))
	bincentres.append((v1+v2)/2.)
	

### above max? ### - put all into maxv bin (centred)
bmax_bool = xyzpd['z'].ge(max(bins))
subset = xyzpd[bmax_bool.values]
sums.append(np.sum(subset['latfac']))
#bincentres.append(maxv+(t/2.))
bincentres.append(maxv)


# integral of cos(x) from -pi/2 to pi/2 = 2.
# integral of y = 1 from -pi/2 to pi/2 = pi
# so the following factorisation should return similar values to giving each count a value of 1 regardless of latitude
sums_orig = np.array(sums)
sums = np.multiply(np.pi/2., sums)

### what histogram should look like
#plt.bar(bincentres, sums)
#plt.show()

tosave = []
for i in range(len(bincentres)):
	tosave.append([bincentres[i], sums[i]])

np.savetxt('myhisto.temp', tosave)




