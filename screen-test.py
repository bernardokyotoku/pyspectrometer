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
device.setRowlandRadius(1E-3) # important to define input position
device.input.setN(1)
device.setInputWidth(3e-6)
device.setInputAngle(180*np.pi/180)
device.setInputCenterOnRowlandCircle()
device.setInputFieldMode(spectrometer.rectModeFunction)
device.setNumberOfGratingGrooves(999)
device.makeConcaveScreen()
#device.viewSpectrometer()
print 'Simulating'
print ' Calculating Eletromagnetic Field on Grating'
device.propagateToGrating(1.5E-6)
device.plotGratingField()
raw_input()
