#!/usr/bin/env python
from numpy import *
from scipy import optimize
import scientific_utils as su
from matplotlib.pyplot import *
import spectrometer as sp
reload(su)
import copy
	
class WavePropagator:
def __init__(spectrometer):
	E=0
	keff=0
	spectrometer=spectrometer
	spectrometer.grating.E = [0]
	setAperturePoints(spectrometer.input,100)
	setCosMode(spectrometer.input)
	spectrometer.input.E = spectrometer.input.mode
	eval_r1=True
	eval_r2=True

def setCosMode(aperture):
	aperture.mode = cos(pi*aperture.x/aperture.width)*sqrt(2/aperture.width)
	
def setAperturePoints(aperture,n):
	aperture.n = n
	spectrometer.setCenter(aperture)
	aperture.eval_r2 = True
	angle = aperture.angle + pi/2
	aperture.x = linspace(-aperture.width/2,aperture.width/2,n)
	aperture.points = transpose(vstack((cos(angle)*aperture.x,sin(angle)*aperture.x)))+aperture.center		
	
def propagateToGrating(keff):
	if keff==keff:
		return
	if eval_r1:
		r1 = spectrometer.grating.points - spectrometer.input.center
		alpha = su.angle(r1)-spectrometer.input.angle
		alphaI = (su.angle(-r1)-spectrometer.grating.angles)
		r1 = sqrt(sum(r1**2,1))
	if len(spectrometer.grating.E) != spectrometer.grating.n:
		spectrometer.grating.E=zeros(spectrometer.grating.n,dtype='complex')
	for i in range(spectrometer.grating.n):
		dr1 = sin(alpha[i])*spectrometer.input.x
		spectrometer.grating.E[i] = sum(spectrometer.input.mode*exp(-1j*keff*dr1))
	dx = spectrometer.input.x[1]-spectrometer.input.x[0]
	spectrometer.grating.E = spectrometer.grating.E*(1+cos(alpha))*exp(-1j*keff*r1)*dx
	spectrometer.grating.E = spectrometer.grating.E/(2*sqrt(2*pi*r1/keff))
	keff = keff
		
def propagateTo(target):
	if not hasattr(target,'E'):
		target.E = zeros(target.n,dtype='complex')
	if len(target.E) != target.n:
		target.E = zeros(target.n,dtype='complex')
	if target.eval_r2:
		target.eval_r2 = False
		target.r2 = zeros((target.n,spectrometer.grating.n))
		target.alphaD = zeros((target.n,spectrometer.grating.n))
		target.sin_diff = zeros((target.n,spectrometer.grating.n))
		target.cos_sum  = zeros((target.n,spectrometer.grating.n))
		for i in range(target.n):
			r  = target.points[i]-spectrometer.grating.points
			target.r2[i] = sqrt(sum(r**2,1))
			target.alphaD[i] = (spectrometer.grating.angles-su.angle(r))
			target.sin_diff[i] = (sin(alphaI)-sin(target.alphaD[i]))/2
			target.cos_sum[i] = (cos(alphaI)+cos(target.alphaD[i]))
	for i in range(target.n):
		EfromEachFacet = sinc(keff*spectrometer.grating.widths*target.sin_diff[i]/pi)
		EfromEachFacet = EfromEachFacet*spectrometer.grating.widths
		EfromEachFacet = EfromEachFacet*exp(-1j*keff*target.r2[i])/sqrt(target.r2[i])
		EfromEachFacet = EfromEachFacet*spectrometer.grating.E*target.cos_sum[i]
		target.E[i] = sum(EfromEachFacet)*sqrt(keff/8/pi)
		
def fractionCoupledInto(waveguide):
	## Calculate the internal product using the simpson rule.
	if waveguide.n%2==0:
		print 'Number of points need to be odd.'
		return -1
	h = waveguide.x[1]-waveguide.x[0]
	weights = [1, 4] + (waveguide.n-3)/2*[2, 4] + [1]
	return (sum(weights*waveguide.mode*abs(waveguide.E))*h/3)**2
	
def transmissionSpectrum(keffs,waveguide):
	spectrum = zeros(len(keffs))
	for i in range(len(keffs)):
		propagateToGrating(keffs[i])
		propagateTo(waveguide)
		spectrum[i] = fractionCoupledInto(waveguide)
	return spectrum

def findPeakTransmissioAngleAt(keff,aperture):
	propagateToGrating(keff)
	N = su.fwhm(abs(spectrometer.grating.E))
	deltaLambdaEff = 2*pi/(spectrometer.order*N*keff)
	a = spectrometer.grating.pitch
	m = spectrometer.order
	lambdaEff = 2*pi/keff
	inputAngle = spectrometer.input.angle
	startAngle = arcsin(m/a*(lambdaEff-deltaLambdaEff/2)+sin(inputAngle))
	endAngle   = arcsin(m/a*(lambdaEff+deltaLambdaEff/2)+sin(inputAngle))
	def f(angle):
		spectrometer.changeApertureAngle(aperture,angle)
		propagateTo(aperture)

		return -fractionCoupledInto(aperture)
	return optimize.golden(f,brack = (startAngle,endAngle), tol=1e-6,full_output=True)
	
def optimizeWidth(aperture):
	propagateToGrating(aperture.keff)
	def f(width):
		aperture.width = width
		setAperturePoints(aperture,21)
		setCosMode(aperture)
		propagateTo(aperture)
		a = fractionCoupledInto(aperture)
		print 'fraction =', a
		return -fractionCoupledInto(aperture)
	return optimize.golden(f,brack = (100e-9,10e-6), tol=1e-6)		

def setOutputArray(keffs):
	aperture = sp.Aperture()
	aperture.width = 1.4e-6
	aperture.angle = 30*pi/180
	spectrometer.setCenter(aperture)
	setAperturePoints(aperture,51)
	setCosMode(aperture)		
	spectrometer.outputs = []
	i=0
	for keff in keffs:
		ape = copy.deepcopy(aperture)
		ape.keff = keff
		findPeakTransmissioAngleAt(keff,ape)
		
		#optimizeWidth(ape)
		spectrometer.outputs.append(ape)
		i+=1
		print i,' done' 
	
	
			
def plotInputField():
	pass

def plotGratingField():
	figure()
	plot(abs(spectrometer.grating.E))
