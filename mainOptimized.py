#! /usr/bin/env python
print 'Loading Modules'
import spectrometer
import RefractiveIndex
import matplotlib.pyplot as plt
import numpy as np
import scientific_utils as su
reload(spectrometer)
print 'Building Spectrometer'
device = spectrometer.Spectrometer()
device.setFilmThickness(250e-9)
device.setRefractiveIndexStackList([RefractiveIndex.nSilica,RefractiveIndex.nSi,RefractiveIndex.nSilica])
device.setPolarization('TE')
device.setRefractiveIndex('ThinFilm')
device.setAberrationFreeWavelength(1550e-9)
outWave = np.linspace(1550e-9,1650e-9,20)#;import pdb;pdb.set_trace()
device.setOutputWavelengths(outWave)
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
device.makeUniformGrating(centralGroove=145)
p = device.grating.points[0]-device.grating.points[-1]
print 'gratingLength=',np.sqrt(sum(p**2))
print 'gratingAngle=',su.angle(p)*180/np.pi
print ' Calculating Grating Blaze Angles'
device.blazeGratingToPoint(device.getAberrationFreePoint())
#fig=figure(1)
print 'Calculating output waveguide position, angle and width'
device.calculateOutputsArray()
device.viewSpectrometer()
print 'Simulating'
print ' Calculating Eletromagnetic Field on Grating'
device.propagateToGrating(device.getAberrationFreeWavelength())
#plt.figure(2)
#device.grating.plotField()
#plt.show()
swavelengths = [np.linspace(w-1E-9,w+1E-9,1000) for w in outWave]
spectra = []
for i in range(len(swavelengths)):
	wavelengths = swavelengths[i]
	output = device.output_list[i]
	spectra.append(device.transmissionSpectrum(wavelengths,output))
	print len(spectra)
	
#figure(3)
#plot(swavelengths,spectrum)
##s.output.points = s.grating.aberrationFreePoint
##s.output.widths = [2e-6]
##s.output = s.makeRowlandCircleEdgePoints(65*pi/180,75*pi/180,3e-6,2e-6)
##s.output = s.makeLineAlongWidth(OutputWidth,optimizedOutputAngle,optimizedPoint,100)
##s.plot()
#
#s.output = spectrometer.Aperture()
#s.output.width = 1.4e-6
#s.output.angle = s.grating.aberrationFreeAngle
#s.setCenter(s.output)
#print ' Calculating Eletromagnetic Output'
#
#s.arc = s.setRowlandCircleEdgePoints(75*pi/180,77*pi/180,2e-6,1e-6)
#wp.propagateTo(s.arc)
#figure(1)
##s.arc.plotField()
#
#wp.setAperturePoints(s.output,51)
##wp.propagateTo(s.output)
##figure(2)
##s.output.plotField()
#
##s.image = spectrometer.Grid(10e-6,10e-6,1e-7,s.output.center)
##s.plot()
##s.image.plot()
##wp.propagateTo(s.image)
##figure(3)
##s.image.plotField()
#
#wp.setCosMode(s.output)
#print 'angle2',s.output.angle
#a= wp.findPeakTransmissioAngleAt(keff(1500e-9),s.output)
#wp.propagateTo(s.output)
#s.output.plotField()
#print 'angle',s.output.angle
##print 'angle2',s.output.angle
#
#
##v = s.grating.points[0:-2]-s.grating.points[1:-1]
##print sqrt(sum(v**2,1))
##print sum(abs(s.grating.E[1:-1])**2*sqrt(sum(v**2,1)) )
#print sum(abs(s.input.mode)**2)*(s.input.x[1]-s.input.x[0])
#print sum(abs(s.output.E)**2)*(s.output.x[1]-s.output.x[0])
#
#s.outputs = s.arc
##s.outputs.points = array([s.output.center])
##s.outputs.angles = array([s.output.angle])
##s.outputs.widths = array([s.output.width])
##s.plot()
#print 'Saving to File,','spectrometer.txt'
#device.exportSpectrometerDetails('spectrometer.txt')


