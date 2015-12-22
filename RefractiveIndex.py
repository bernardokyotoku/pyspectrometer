#!/usr/bin/env python
import numpy

#Silicon
def nSi(wavelength):
	"""
	Return the refractive index for the given wavelength in meters.

	Source: Tropf,Harris,Thomas -  Optical materials visible and infrared  
	"""
	l = [0.301516485**2,1.13475115**2,1104.0**2]
	A = [10.6684293,0.00304347484,1.54133408]
	w = (wavelength*1e6)**2
	return numpy.sqrt(1+(A[0]/(w-l[0])+A[1]/(w-l[1])+A[2]/(w-l[2]))*w)

#Fused Silica
def nSilica(wavelength):
	"""
	Return the refractive index for the given wavelength in meters.

	Source: Tropf,Harris,Thomas -  Optical materials visible and infrared  
	"""
	l = [0.0684043**2,0.1162414**2,9.896161**2]
	A = [0.6961663,0.4079426,0.8974794]
	w = (wavelength*1e6)**2
	return numpy.sqrt(1+(A[0]/(w-l[0])+A[1]/(w-l[1])+A[2]/(w-l[2]))*w)
  
#Silicon Nitride
def nSi3N4(wavelength):
        """ 
        Return the refractive index for the given wavelength in meters.

	Source: Someone from the cornell nanophotonics group
        """
        l = [19885.02553,-16180.13915]
        A = [0.71839,0.2749]
        w = (wavelength*1e9)**2
        return 1+(A[0]/(w-l[0])+A[1]/(w-l[1]))*w

