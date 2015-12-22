#!/usr/bin/env python
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def LD (R,m,a,wavelength,angleIn): 
	return R*m/(a*(1+(m*wavelength/a+np.sin(angleIn))**2))

R = 750
m = -10
a = 4
wavelength = 1.5

N = 100
angles = np.linspace(0,np.pi/2,N)
#import pdb;pdb.set_trace()
plt.plot(angles,np.abs(LD(R,m,a,wavelength,angles)))
plt.ylabel('Linear Dispersion $\\mu m/\\mu m$')
plt.xlabel('$\\theta_{in}$ degrees')
plt.xticks(np.arange(0,np.pi/2,np.pi/12),[str(i) for i in range(0,90,15)])
plt.show()
