import numpy as np
from scipy import optimize
import ConfigParser
import scientific_utils as su
import copy
try:
	import matplotlib.pyplot as plt
except ImportError:
	print 'matplotlib not found, visualization methods disabled'
	plt = 0

class Aperture:
	def __init__(self):
		self.eval_r2 = True
		self.n = 21
		self.points=0
		self.widths=0
		self.angles=0
	
	def setAngle(self,angle):
		self.angle = angle

	def setWidth(self,width):
		self.width = width

	def setCenter(self,center):
		self.center = center

	def setN(self,N):
		self.n = N
		
	def makePoints(self):
		self.x = np.linspace(-self.width/2,self.width/2,self.n)
		self.points = (su.unitVector(self.angle+np.pi/2).transpose()*self.x).transpose() + self.center

	def pointsParametersOK(self):
		return dir(self).__contains__('width') and dir(self).__contains__('n') and dir(self).__contains__('angle') and 	dir(self).__contains__('center')
		
	def plot(self):
		if plt==0:
			print 'matplotlib not loaded.'
			return -1
		
		if np.size(self.points)==0:
			return
		plt.plot(self.points[...,0],self.points[...,1],'o')
		vy = np.array([np.cos(self.angle),np.sin(self.angle)])
		vx = np.array([np.cos(self.angle-np.pi/2),np.sin(self.angle-np.pi/2)])
		arrow = np.zeros((8,2))
		arrow[0] = self.center+vx*0.25*self.width
		arrow[1] = self.center-vx*0.25*self.width
		arrow[2] = self.center-vx*0.25*self.width+vy*0.5*self.width
		arrow[3] = self.center-vx*0.50*self.width+vy*0.5*self.width
		arrow[4] = self.center+vy*1*self.width
		arrow[5] = self.center+vx*0.50*self.width+vy*0.5*self.width
		arrow[6] = self.center+vx*0.25*self.width+vy*0.5*self.width
		arrow[7] = arrow[0]
		plt.plot(arrow[...,0],arrow[...,1])

	def plotField(self):
		if plt==0:
			print 'matplotlib not loaded.'
			return -1
		
		if np.size(self.points)==0:
			return
		plt.plot(self.x,np.abs(self.E))

class DiffractionGrating:
	def __init__(self):
		pass
	
	def plot(self):
		if plt==0:
			print 'matplotlib not loaded.'
			return -1
		plt.plot(self.points[:,0],self.points[:,1],'o')
		
	def plotField(self):
		if plt==0:
			print 'matplotlib not loaded.'
			return -1
		plt.plot(np.abs(self.E))		
	
class Input:
	def __init__(self,width):
		self.changed = True
		self.width = width
		self.point=0
		self.angle=0

	def plot(self):
		if plt==0:
			print 'matplotlib no loaded.'
			return -1
		plt.plot([self.point[0]],[self.point[1]],'o')
	

class Grid:
	def __init__(self,Xlength,Ylength,pitch,center):
		X = np.arange(-Xlength/2,Xlength/2,pitch)
		Y = np.arange(-Ylength/2,Ylength/2,pitch)
		self.eval_r2 = True
		self.nX = len(X)
		self.nY = len(Y)
		self.n = self.nY*self.nX
		m = np.array(meshgrid(X,Y))
		self.points =  m.flatten().reshape(2,self.n).transpose()+center
	
	def plotField(self):
		if plt==0:
			print 'matplotlib not loaded.'
			return -1

		plt.imshow(flipud(np.abs(self.E).reshape(self.nY,self.nX)))

	def plot(self):
		if plt==0:
			print 'matplotlib not loaded.'
			return -1
		
		if size(self.points)==0:
			return
		plt.plot(self.points[...,0],self.points[...,1],'o')
class Spectrometer:
	def __init__(self):
		self.input = Aperture()
		self.output = Aperture()
		self.grating = DiffractionGrating()
		self.output_list = []
		self.outputs = Aperture()
		self.keff = 0
		self.eval_r1 = True
		self.eval_r2 = True
 	
	def setThinFilmNeff(self,mode):
		if dir(self).__contains__('thickness') and dir(self).__contains__('material_stack') and dir(self).__contains__('polarization'):
			import thinfilms
			self.neff = lambda w: thinfilms.neff(w,self.thickness,self.material_stack[0],self.material_stack[2],self.material_stack[1],self.polarization)[mode-1]
	
	def setFilmThickness(self,Thickness):
		self.thickness = Thickness
	
	def setRefractiveIndexStackList(self,List):
		self.material_stack = List
	
	def setPolarization(self,pol):
		self.polarization = pol
	
	def setRowlandRadius(self,RowlandRadius):
		self.radius = RowlandRadius
	
	def setDiffractionOrder(self,Order):
		self.order = Order
	
	def setGratingPitch(self,pitch):
		self.grating.pitch = pitch
	
	def setGratingFacetWidth(self,width):
		self.grating.width = width
	
	def setNumberOfGratingGrooves(self,n):
		self.grating.n = n
	
	def setAberrationFreeWavelength(self,wavelength):
		self.grating.AberrationFreeWavelength = wavelength
		self.grating.AberrationFreeKeff = self.convertToKeff(wavelength)
	
	def setRefractiveIndex(self,n,mode=1):
		if type(n).__name__== 'str':
			if n == 'ThinFilm':
				self.setThinFilmNeff(mode)
		elif type(n).__name__ == 'function':
			self.neff = n
		elif type(n).__name__ == 'int':
			self.neff = lambda w : float(n)*w
		elif type(n).__name__ == 'float':
			self.neff = lambda w : n*w
	
	def convertToKeff(self,wavelengths):
		def isNumber(num):
			numberTypes = np.array(['int','float','long','complex','float64'])
			table = type(num).__name__ == numberTypes
			return sum(table)>0
		if isNumber(wavelengths):
			return 2*np.pi*self.neff(wavelengths)/wavelengths
		else:
			return np.array([ 2*np.pi*self.neff(wavelength)/wavelength for wavelength in wavelengths])
	
	def setApertureCenterOnRowlandCircle(self,aperture):
		aperture.setCenter(self.radius*su.unitVector(2*aperture.angle-2*np.pi))

	def setOutputWavelengths(self,wavelengths):
		self.outputs.keffs = self.convertToKeff(wavelengths) 
		self.outputs.wavelengths = wavelengths
	
	def setOutputWidth(self,width):
		self.output.setWidth(width)
		if self.output.pointsParametersOK():
			self.output.makePoints()

	def setOutputAngle(self,angle):
		self.output.setAngle(angle)
		if self.output.pointsParametersOK():
			self.output.makePoints()

	def setOutputFieldMode(self,modeFunction):
		self.output.mode = modeFunction
	
	def setOutputCenterOnRowlandCircle(self):
		self.output.setCenter(self.radius*su.unitVector(2*self.output.angle-2*np.pi))
		if self.output.pointsParametersOK():
			self.output.makePoints()

	def setInputWidth(self,width):
		self.input.setWidth(width)
		if self.input.pointsParametersOK():
			self.input.makePoints()
	
	def setInputAngle(self,angle):
		self.input.setAngle(angle)
		if self.input.pointsParametersOK():
			self.input.makePoints()
	
	def setInputFieldMode(self,modeFunction):
		self.input.mode = modeFunction

	def setInputCenterOnRowlandCircle(self):
		self.input.setCenter(self.radius*su.unitVector(2*self.input.angle-2*np.pi))
		if self.input.pointsParametersOK():
			self.input.makePoints()
	
	def setRowlandCircleEdgePoints(self,startAngle,endAngle,stepSize,width):
		output = Aperture()
		output.eval_r2 = True
		output.angles = np.linspace(startAngle,endAngle,self.radius*(endAngle-startAngle)/stepSize)
		output.n = len(output.angles)
		output.points = np.array([np.cos(output.angles),np.sin(output.angles)]).transpose()*self.radius
		output.x = output.angles*180/np.pi
		output.widths = ones(output.n)*width
		output.angles = output.angles/2+np.pi
		return output

	def getAberrationFreeWavelength(self):
		return self.grating.AberrationFreeWavelength	
	
	def getAberrationFreePoint(self):
		self.grating.AberrationFreePoint = np.array([np.cos(2*self.grating.AberrationFreeAngle),np.sin(2*self.grating.AberrationFreeAngle)])*self.radius
		return self.grating.AberrationFreePoint
	
	def makeUniformGrating(self,centralGroove=None):
		self.changed = True
		self.grating.AberrationFreeAngle = np.arcsin(np.sin(self.input.angle)+self.order/self.grating.pitch*2*np.pi/self.grating.AberrationFreeKeff)
		self.grating.AberrationFreePoint = np.array([np.cos(2*self.grating.AberrationFreeAngle),np.sin(2*self.grating.AberrationFreeAngle)])*self.radius
		self.grating.AberrationFreeAngle += np.pi
		if centralGroove == None:
			centralGroove = self.grating.n/2	
		Y = np.arange(-self.grating.n/2,self.grating.n/2,1)*self.grating.pitch
		print "len Y ",len(Y)
		print "graing.n ",self.grating.n
		angle = np.arcsin(Y/self.radius)
		X = self.radius*(1-2*np.cos(angle))
		self.grating.points = np.vstack((X,Y)).T
	
	def makeOneStigmaticPointGrating(self,centralGroove=None):
		self.changed = True
		self.grating.AberrationFreeAngle = np.arcsin(np.sin(self.input.angle)+self.order/self.grating.pitch*2*np.pi/self.grating.AberrationFreeKeff)
		self.grating.AberrationFreePoint = np.array([np.cos(2*self.grating.AberrationFreeAngle),np.sin(2*self.grating.AberrationFreeAngle)])*self.radius
		self.grating.AberrationFreeAngle += np.pi
		if centralGroove == None:
			centralGroove = self.grating.n/2	
		def normalizedPathLength(theta):
			I = self.input.center[0]
			S = self.grating.AberrationFreePoint
			Gx = (1-2*np.cos(theta))*self.radius 
			Gy = 2*np.sin(theta)*self.radius
			IG = np.sqrt((Gx-I[0])**2+(Gy-I[1])**2)	
			GS = np.sqrt((Gx-S[0])**2+(Gy-S[1])**2)
			return self.grating.AberrationFreeKeff*(IG+GS)/2/np.pi
		Ref = normalizedPathLength(0)
		self.grating.points = np.zeros((self.grating.n,2))
		for i in np.arange(self.grating.n):
			def g(theta):
				s = normalizedPathLength(theta)
				s = s-self.order*(i-centralGroove)
				s = s-Ref
				return s
			angle = optimize.fsolve(g,0)
			self.grating.points[i,:] = np.array([(-2*np.cos(angle)+1),2*np.sin(angle)])*self.radius	
	
	def blazeGratingToPoint(self,blazePoint):
		changed = True
		I = self.input.center-self.grating.points
		mI = np.sqrt(np.sum(I**2,1))
		S = blazePoint-self.grating.points
		angleI = su.angle(I)
		self.grating.angles = (angleI+su.angle(S))/2
		an = np.diff((self.grating.angles[1:]+self.grating.angles[:-1])/2)
		an = np.hstack((an,an[-1]))
		an = np.hstack((an[0],an))
		self.grating.widths = an*mI/np.cos(self.grating.angles-angleI)
	
	def viewSpectrometer(self):
		if plt==0:
			print 'matplotlib not loaded'
			return -1
		self.grating.plot()
		self.input.plot()
		if hasattr(self,'output_list'):
			for output in self.output_list:
				output.plot()
		t=np.linspace(0,2*np.pi,100)
		plt.plot(self.radius*np.cos(t),self.radius*np.sin(t))
		plt.axes().set_aspect('equal')
		plt.show()
		
	def write(self,filePointer):
		cfg = ConfigParser.ConfigParser()
		cfg.add_section('Parameters')
		cfg.set('Parameters','radius',self.radius)
		cfg.set('Parameters','order',self.order)
		cfg.set('Parameters','aberrationFreeKeff',self.grating.AberrationFreeKeff)
		cfg.set('Parameters','aberrationFreeWavelength',self.grating.AberrationFreeWavelength)
		cfg.add_section('Detail')
		cfg.set('Detail','GratingPoints',repr(self.grating.points)	)
		cfg.set('Detail','GratingAngles',repr(self.grating.angles.view()	))
		cfg.set('Detail','GratingWidths',repr(self.grating.widths.view()	))
		cfg.set('Detail','InputPoint'	,repr(self.input.center.view()		))
		cfg.set('Detail','InputAngle'	,repr(self.input.angle				))
		cfg.set('Detail','InputWidth'	,repr(self.input.width				))
		points = np.array([output.center[0] for output in self.output_list])
		angles = np.array([output.angle for output in self.output_list])
		widths = np.array([output.width for output in self.output_list])
		cfg.set('Detail','OutputPoints' ,repr(points.view()		))
		cfg.set('Detail','OutputAngles'	,repr(angles.view()		))
		cfg.set('Detail','OutputWidths'	,repr(widths.view()		))
		cfg.write(filePointer)
		
	def propagateToGrating(self,wavelength):
		if self.keff== self.convertToKeff(wavelength):
			return
		self.keff = self.convertToKeff(wavelength)
		if self.eval_r1:
			self.eval_r1 = False
			self.r1 = (self.grating.points - self.input.center).T
			ar1 = np.arctan2(self.r1[1],self.r1[0])
			alpha = -self.input.angle
			print 'angles',self.grating.angles.shape
			print 'ar1',ar1.shape
			print 'r1',self.r1.shape
			self.alphaI = (ar1+np.pi-self.grating.angles)
			self.r1 = np.sqrt(np.sum(self.r1**2,0))
			print 'sr1',self.r1.shape
			self.grating.E=np.zeros(self.grating.n,dtype='complex')
			self.input.E = self.input.mode(self.input.x,self.input.width)
		x = self.input.x
		E = self.input.E
#		for i in range(self.grating.n):
#			dr1 = np.sin(alpha[i])*self.input.x
#			self.grating.E[i] = np.sum(self.input.E*np.exp(-1j*self.keff*dr1))
		if len(E.shape) < 2:
			E.shape += (1,)
		if len(x.shape) < 2:
			x.shape += (1,)
		U = self.grating.E
		k = self.keff
		r1 = self.r1
		dx = (self.input.x[1]-self.input.x[0]).reshape(())
		U = E*np.exp(-1j*k*np.sin(alpha)*x)
		U = np.sum(U,0)
		U = U*(1+np.cos(alpha))*np.exp(-1j*k*r1)*dx/(2*np.sqrt(2*np.pi*r1/k))
		print 'U',U.shape
		self.grating.E = U
			
	def propagateTo(self,target):
		if not hasattr(target,'E'):
			target.E = np.zeros(target.n,dtype='complex')
		if len(target.E) != target.n:
			target.E = np.zeros(target.n,dtype='complex')
		if target.eval_r2:
			target.eval_r2 = False
			target.points.shape = (1,)+ target.points.shape
			self.grating.points.shape += (1,)
			O = target.points.transpose((2,0,1))
			G = np.transpose(self.grating.points,(1,0,2))
			R  = O-G
			target.r2 = np.sqrt(np.sum(R**2,0)).T 
			angles = self.grating.angles
			target.alphaD = angles-np.arctan2(R[1],R[0]).T
			target.sin_diff = (np.sin(self.alphaI)-np.sin(target.alphaD))/2
			target.cos_sum = (np.cos(self.alphaI)+np.cos(target.alphaD))
#			print 'target',O.shape
#			print 'grating',G.shape
#			print 'r',r.shape
#			print 'R',R.shape
#			print 'angle',self.grating.angles.shape
#			print 'alphaD',alphaD.shape
#			print 'alphaI',self.alphaI.shape
#			print 'sin_diff',sin_diff.shape
#			print 'cos_sum',cos_sum.shape	
#			for i in range(target.n):
#				target.r2[i] = np.sqrt(np.sum(r**2,1))
#				target.alphaD[i] = (self.grating.angles-su.angle(r))
#				target.sin_diff[i] = (np.sin(self.alphaI)-np.sin(target.alphaD[i]))/2
#				target.cos_sum[i] = (np.cos(self.alphaI)+np.cos(target.alphaD[i]))
		for i in range(target.n):
			r2 = target.r2
			sin_diff = target.sin_diff
			cos_sum = target.cos_sum
			alphaD = target.alphaD
			k = self.keff
			w = self.grating.widths
			E = self.grating.E
			print 'w',w.shape
			U = np.sum(w*np.sinc(k*w*sin_diff/np.pi)*np.exp(-1j*k*r2)/np.sqrt(r2)*E*cos_sum,0)
			print 'U',U.shape
			EfromEachFacet = np.sinc(self.keff*self.grating.widths*target.sin_diff[i]/np.pi)
			EfromEachFacet = EfromEachFacet*self.grating.widths
			EfromEachFacet = EfromEachFacet*np.exp(-1j*self.keff*target.r2[i])/np.sqrt(target.r2[i])
			EfromEachFacet = EfromEachFacet*self.grating.E*target.cos_sum[i]
			target.E[i] = np.sum(EfromEachFacet)*np.sqrt(self.keff/8/np.pi)
			

	def fractionCoupledInto(self,aperture):
		## Calculate the internal product using the simpson rule.
		if aperture.n%2==0:
			print 'Number of points need to be odd.'
			return -1
		h = aperture.x[1]-aperture.x[0]
		weights = [1, 4] + (aperture.n-3)/2*[2, 4] + [1]
		return (np.sum(weights*aperture.mode(aperture.x,aperture.width)*np.abs(aperture.E))*h/3)**2
		
	def transmissionSpectrum(self,wavelengths,waveguide):
		spectrum = np.zeros(len(wavelengths))
		for i in range(len(wavelengths)):
			self.propagateToGrating(wavelengths[i])
			self.propagateTo(waveguide)
			spectrum[i] = self.fractionCoupledInto(waveguide)
		return spectrum
	
	def findPeakTransmissionAngleAt(self,wavelength,aperture):
		self.propagateToGrating(wavelength)
		print 'Ed',np.abs(self.grating.E)
		N = su.fwhm(np.abs(self.grating.E))
		deltaLambdaEff = wavelength/(self.order*N)/self.neff(wavelength)
		a = self.grating.pitch
		m = self.order
		lambdaEff = wavelength/self.neff(wavelength)
		inputAngle = self.input.angle
		startAngle = np.arcsin(m/a*(lambdaEff-deltaLambdaEff/2)+np.sin(inputAngle))
		endAngle   = np.arcsin(m/a*(lambdaEff+deltaLambdaEff/2)+np.sin(inputAngle))
		def f(angle):
			aperture.setAngle(angle+np.pi)
			self.setApertureCenterOnRowlandCircle(aperture)
			aperture.makePoints()
			self.propagateTo(aperture)
			return -self.fractionCoupledInto(aperture)
		return optimize.golden(f,brack = (startAngle,endAngle), tol=1e-6,full_output=True)
		
	def optimizeWidth(self,aperture):
		self.propagateToGrating(aperture.wavelength)
		def f(width):
			aperture.width = width
			setAperturePoints(aperture,21)
			setCosMode(aperture)
			propagateTo(aperture)
			a = fractionCoupledInto(aperture)
			print 'fraction =', a
			return -fractionCoupledInto(aperture)
		return optimize.golden(f,brack = (100e-9,10e-6), tol=1e-6)		
	
	def calculateOutputsArray(self):
		for i in range(len(self.outputs.wavelengths)):
			aperture = copy.deepcopy(self.output)
			aperture.keff = self.outputs.keffs[i]
			aperture.wavelength = self.outputs.wavelengths[i]
			self.findPeakTransmissionAngleAt(aperture.wavelength,aperture)
			self.output_list.append(aperture)
			print i,' done' 
		
	def plotInputField(self):
		pass
	
	def plotGratingField(self):
		figure()
		plot(np.abs(self.grating.E))
		
	def exportSpectrometerDetails(self,filename):
		fd = open(filename,'w')
		self.write(fd)
		fd.close()


def cosModeFunction(x,width):
	return np.cos(np.pi*x/width)*np.sqrt(2/width)
