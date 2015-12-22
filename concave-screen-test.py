#! /usr/bin/env python
print 'Loading Modules'
import spectrometer
import RefractiveIndex
import matplotlib.pyplot as plt
import numpy as np
import scientific_utils as su
reload(spectrometer)
plt.rc('text',usetex=True)
plt.rc('font', family='Serif')
plt.rc('legend',fontsize=10)
print 'Building Spectrometer'
toDegree = lambda x:x*180/pi
c = 'red'
for n in [21,101]:
	device = spectrometer.Spectrometer()
	device.setFilmThickness(250e-9)
	device.setRefractiveIndexStackList([RefractiveIndex.nSilica,RefractiveIndex.nSi,RefractiveIndex.nSilica])
	device.setPolarization('TE')
	device.setRefractiveIndex('ThinFilm')
	device.setRowlandRadius(1E-3) # important to define input position
	device.input.setN(n)
	device.setInputWidth(3e-6)
	device.setInputAngle(180*np.pi/180)
	device.setInputCenterOnRowlandCircle()
	device.setInputFieldMode(spectrometer.rectModeFunction)
	device.setNumberOfGratingGrooves(50)
	device.makeConcaveScreen()
	#device.viewSpectrometer()
	print 'Simulating'
	print ' Calculating Eletromagnetic Field on Grating'
	device.propagateToGrating(1.5E-6)
	pi = np.pi
	p = -(device.grating.points-device.input.center)
	angle = lambda x:np.arctan2(x[:,1],x[:,0])
	alpha = angle(p)
	plt.scatter(toDegree(alpha),abs(device.grating.E),color=c,label='N %d'%n)
	c = 'green'
#device.plotGratingField()
sin = np.sin
exp = np.exp
cos = np.cos
sinc = np.sinc
sqrt = np.sqrt
#import pdb;pdb.set_trace()
k = device.keff
r = sqrt(np.sum(p**2,1))
w = device.input.width
arg = 0.5*k*sin(alpha)*w
E = 0.5*sqrt(k/2/pi)*(1+cos(alpha))*w*sinc(arg/pi)/sqrt(r)
O0 = (1+cos(alpha))*w*sinc(arg/pi)
O1 = (sinc(arg/pi)-cos(arg))*w*cos(alpha)*3/(2*r*k)
E = 0.5*sqrt(k/2/pi)/sqrt(r)*(O0+O1)
plt.plot(toDegree(alpha),abs(E),label="Analytic")
plt.gcf().set_size_inches([3.2,2.4])
plt.subplots_adjust(bottom=0.2,left=0.20,right=0.95,top=0.95)
plt.legend(loc='lower left')
plt.ylim([0,0.1])
plt.xlim([-15,15])
plt.xlabel('Angle (degrees)')
plt.ylabel('Field amplitude')
plt.savefig('compare-analytic.pdf')
print 'k=%e'%k
print 'R=%e'%np.mean(r)
raw_input()
