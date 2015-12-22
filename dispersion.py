#!/usr/bin/env python

from thinfilms import neff
from RefractiveIndex import nSi,nSilica
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


N=100
n = lambda w,t: neff(w,t*1E-9,nSilica(w),nSilica(w),nSi(w),pol='TE')

wavelengths = np.linspace(1450E-9,1600E-9,N)
wavem = (wavelengths[1:]+wavelengths[:-1])/2
neffs = []
for thickness in range(200,301,50):
	nf = np.array([n(w,thickness)[0] for w in wavelengths]).T
	neffs.append(nf)
	ng = (nf[1:]+nf[:-1])/2 - wavem*np.diff(nf)/np.diff(wavelengths)
	plt.plot(wavem*1E9,ng,label=str(thickness)+' nm')
plt.subplots_adjust(right=0.94)
plt.xlabel('Wavelength (nm)')
plt.ylabel('$n_g$')
plt.legend(loc='lower left',title='Si film thickness')
plt.xlim([1450,1600])
ran = range(1450,1601,25)
plt.xticks(ran,[str(i) for i in ran])
plt.text(1413,3.775,'b)')
plt.gcf().set_size_inches([3.4,2.4])
plt.savefig('slabGroupDispersion.pdf')
#plt.show()

fig2 = plt.figure(2)
i=0
for nf in neffs:
	thick = range(200,301,50)[i]
	plt.plot(wavelengths*1E9,nf,label=str(thick)+' nm')
	i+=1
plt.subplots_adjust(right=0.94)
plt.legend(title='Si film thickness')
plt.xlim([1450,1600])
plt.xlabel('Wavelength (nm)')
plt.ylabel('$n_{eff}$')
plt.text(1413,3.09,'a)')
ran = range(1450,1601,25)
plt.xticks(ran,[str(i) for i in ran])
fig2.set_size_inches([3.4,2.4])
plt.legend(loc='lower left',title='Si film thickness')
plt.savefig('slabDispersion.pdf')
#plt.show()
