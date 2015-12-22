#!/usr/bin/env/ python
from numpy import *
from scipy import optimize

def TEEigenValueEquation(neff,nc,ns,nf,T,k0):
	"""
	Try guess ns<neff<nf
	"""
	gc=sqrt(neff**2-nc**2)
	gs=sqrt(neff**2-ns**2)
	kf=sqrt(nf**2-neff**2)
	a=kf*k0*T
	
	return (kf**2-gc*gs)*sin(a)-kf*(gs+gc)*cos(a)

def TMEigenValueEquation(neff,nc,ns,nf,T,k0):
        """
        Try guess ns<neff<nf
        """
        gc=sqrt(neff**2-nc**2)
        gs=sqrt(neff**2-ns**2)
        kf=sqrt(nf**2-neff**2)
        a=kf*k0*T
        
        return (nc**2*ns**2*kf**2-nf**4*gc*gs)*sin(a)-nf**2*kf*(nc**2*gs+ns**2*gc)*cos(a)

def neff(wavelength,thickness,nc,ns,nf,pol='TE'):
	"""
	Calculate the effective refrative index of a film.

	Use ns<nc<nf
	"""
	if callable(nc):
		nc=nc(wavelength)
	if callable(ns):
		ns=ns(wavelength)
	if callable(nf):
		nf=nf(wavelength)
	k0=2*pi/wavelength
	if pol=='TE':
		f = lambda x: TEEigenValueEquation(x,nc,ns,nf,thickness,k0)
	else:
		f = lambda x: TMEigenValueEquation(x,nc,ns,nf,thickness,k0)
	kfMax=sqrt(nf**2-ns**2)
	aMax=k0*thickness*kfMax
	n=4*aMax/pi
	AS = linspace(1e-4,aMax*0.999999,n)
	neffs = sqrt(nf**2-(AS/(thickness*k0))**2)
	neff=[]
	for i in range(len(neffs)-1):
		if sign(f(neffs[i]))!=sign(f(neffs[i+1])):
			n=optimize.brentq(f,neffs[i],neffs[i+1])
			neff.append(n)
	return neff
