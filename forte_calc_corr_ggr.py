
import pyshtools as pysh
import numpy as np

import spectrum
import cross_spectrum

#from pyshtools.spectralanalysis import SHAdmitCorr

### CHECK NORMALIZATION? SEE https://shtools.github.io/SHTOOLS/python-shcoeffs.html
### https://shtools.github.io/SHTOOLS/python-shcoeffs.html

try:		
    flm1 = pysh.SHCoeffs.from_file('t1.temp', format='dov')
    flm2 = pysh.SHCoeffs.from_file('t2.temp', format='dov')
except:
    print('dov format not recognised, using shtools format')
    flm1 = pysh.SHCoeffs.from_file('t1.temp', format='shtools')
    flm2 = pysh.SHCoeffs.from_file('t2.temp', format='shtools')

power1 = flm1.spectrum() #- tot per l, not avg per lm - see , unit=
power2 = flm2.spectrum() #- tot per l, not avg per lm - see , unit=
# https://shtools.github.io/SHTOOLS/spectrum.html

# to get lon/lat, use e.g.
#grid = flm1.expand()
#grid.plot()

lmax = flm1.lmax

admit1, corr1 = flm1.admitcorr(flm2, lmax=lmax)
admit2, corr2 = flm2.admitcorr(flm1, lmax=lmax)

a1 = [ i[0] for i in admit1 ]
err1 = [ i[1] for i in admit1 ]

a2 = [ i[0] for i in admit2 ]
err2 = [ i[1] for i in admit2 ]

#flm1.plot_admitcorr(flm2, lmax=lmax)
#flm2.plot_admitcorr(flm1, lmax=lmax)

#print(corr1)

# calculate correlation coefficient for entire spectrum, r_l
sgg = spectrum.spectrum(flm1.coeffs,lmax=lmax)
shh = spectrum.spectrum(flm2.coeffs, lmax=lmax)
sgh = cross_spectrum.cross_spectrum(flm1.coeffs,flm2.coeffs,lmax=lmax)
rltot = sum(sgh)/(np.sqrt(sum(sgg)) * np.sqrt(sum(shh)))

sgg10 = spectrum.spectrum(flm1.coeffs,lmax=10)
shh10 = spectrum.spectrum(flm2.coeffs, lmax=10)
sgh10 = cross_spectrum.cross_spectrum(flm1.coeffs,flm2.coeffs,lmax=10)
rl10 = sum(sgh10)/(np.sqrt(sum(sgg10)) * np.sqrt(sum(shh10)))

#print(rltot, rl10)

np.savetxt('admit1.temp', a1)
np.savetxt('err1.temp', err1)
np.savetxt('corr1.temp', corr1)
np.savetxt('power1.temp', power1)

np.savetxt('admit2.temp', a2)
np.savetxt('err2.temp', err2)
np.savetxt('corr2.temp', corr2)
np.savetxt('power2.temp', power2)

f = open("rltot.temp", "w")
f.write(str(rltot))
f.close()

f = open("rl10.temp", "w")
f.write(str(rl10))
f.close()

