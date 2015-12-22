#! /usr/bin/env python
print 'Loading Modules'
import spectrometer
import RefractiveIndex
import matplotlib.pyplot as plt
import numpy as np
import scientific_utils as su
reload(spectrometer)
print 'Building Spectrometer'
plt.figure(1, figsize=(3.2,2.4))
plt.rc('text',usetex=True)
plt.rc('font',family='Serif')
plt.rc('legend',fontsize=10)
ax = plt.axes([0.18, 0.16, 0.80, 0.80])
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission (dB)')
for amplitude in [30E-9,100E-9,200E-9]:
	device = spectrometer.Spectrometer()
	device.setFilmThickness(250e-9)
	device.setRefractiveIndexStackList([RefractiveIndex.nSilica,RefractiveIndex.nSi,RefractiveIndex.nSilica])
	device.setPolarization('TE')
	device.setRefractiveIndex('ThinFilm')
	device.setAberrationFreeWavelength(1550e-9)
	device.setOutputWavelengths(np.linspace(1550e-9,1550e-9,1))
	device.setOutputWidth(1.4e-6)
	device.setOutputAngle((30+180)*np.pi/180)
	device.setOutputFieldMode(spectrometer.cosModeFunction)
	device.setRowlandRadius(750e-6)
	device.setDiffractionOrder(10)
	device.setInputWidth(1.2e-6)
	device.setInputAngle((44+180)*np.pi/180)
	device.setInputCenterOnRowlandCircle()
	device.setInputFieldMode(spectrometer.cosModeFunction)
	device.setNumberOfGratingGrooves(340)
	device.setGratingPitch(4e-6)
	print ' Building Aberration Reduced Grating'
	device.makeOneStigmaticPointGrating(centralGroove=145)
	p = device.grating.points[0]-device.grating.points[-1]
	print 'gratingLength=',np.sqrt(sum(p**2))
	print 'gratingAngle=',su.angle(p)*180/np.pi
	print ' Calculating Grating Blaze Angles'
	device.blazeGratingToPoint(device.getAberrationFreePoint())
	device.grating.offLineGrooves(amplitude)
	#fig=figure(1)
	print 'Calculating output waveguide position, angle and width'
	device.calculateOutputsArray()
	#device.viewSpectrometer()
	print 'Simulating'
	print ' Calculating Eletromagnetic Field on Grating'
	device.propagateToGrating(device.getAberrationFreeWavelength())
	#plt.figure(2)
	#device.grating.plotField()
	#plt.show()
	swavelengths = np.linspace(1525e-9,1575e-9,1000)
	output = device.output_list[0]
	spectrum =device.transmissionSpectrum(swavelengths,output) 
	plt.plot(swavelengths*1E9,10*np.log10(spectrum),label='%d nm'%(np.round(amplitude*1E9)))
	plt.xlim(1525,1575)
plt.ylim(-40,0)
plt.legend(title='Amplitude')
plt.savefig('random-offline-test.pdf')


