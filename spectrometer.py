import numpy as np
from scipy import optimize
import ConfigParser
import scientific_utils as su
import copy
try:
	import matplotlib.pyplot as plt
	plt.ion()
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
		'''
		Set the angle of the normal of the aperture. With the normal pointing twards the grating.
		'''
		self.angle = angle
		self.eval_r2 = True

	def setWidth(self,width):
		self.width = width

	def setCenter(self,center):
		self.center = center
		self.eval_r2 = True

	def setN(self,N):
		self.n = N
		
	def makePoints(self):
		self.x = np.linspace(-self.width/2,self.width/2,self.n)#.reshape((self.n,1))
		self.points = (su.unitVector(self.angle+np.pi/2).transpose()*self.x).transpose() + self.center
		self.eval_r2 = True

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

	def offLineGrooves(self,amplitude,distribuition=np.random.rand):
		displacements = (distribuition(self.n)-0.5)*amplitude
		displacements.shape += (1,)
		groove_normals = np.vstack((np.cos(self.angles),np.sin(self.angles))).T
		#import pdb;pdb.set_trace()
		self.points += groove_normals*displacements
	
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
		m = np.array(np.meshgrid(X,Y))
		self.points =  m.flatten().reshape(2,self.n).transpose()+center
	
	def plotField(self):
		if plt==0:
			print 'matplotlib not loaded.'
			return -1

		plt.imshow(np.flipud(np.abs(self.E).reshape(self.nY,self.nX)))

	def plot(self):
		if plt==0:
			print 'matplotlib not loaded.'
			return -1
		
		if len(self.points)==0:
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
		self.input.E = self.input.mode(self.input.x,self.input.width)

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
		angle = np.arcsin(0.5*Y/self.radius)
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
			angle = optimize.fsolve(g,0)[0]
			a = np.array([(-2*np.cos(angle)+1),2*np.sin(angle)])*self.radius	
			#from IPython import embed;embed()
			self.grating.points[i] = a

	def makeConcaveScreen(self,width=1E-3):
		"""
		width - screen width
		posistion - screen center position
		"""
		R = self.radius
		self.changed = True
		N = self.grating.n
		dw = width/N 
		i = np.linspace(-(N-1)/2,(N-1)/2,N)
		y = i*dw
		x = -np.sqrt((2*R)**2-y**2)
		#import pdb;pdb.set_trace()
		self.grating.points = np.vstack((x+R,y)).T
		self.grating.angles = np.arctan2(y,x)
		self.grating.widths = dw*np.ones(N)/np.cos(self.grating.angles)

	def makeScreen(self,width=1E-3,position=(-1E-3,0)):
		"""
		width - screen width
		posistion - screen center position
		"""
		self.changed = True
		N = self.grating.n
		dw = width/N 
		i = np.linspace(-(N-1)/2,(N-1)/2,N)
		i.shape += (1,)
		self.grating.points = position + i*(0,1)*dw	
		self.grating.widths = dw*np.ones(N)
		self.grating.angles = np.zeros(N)

	
	def blazeGratingToPoint(self,blazePoint):
		changed = True
		I = self.input.center-self.grating.points
		mI = np.sqrt(np.sum(I**2,1))
		S = blazePoint-self.grating.points
		angleI = np.arctan2(I[:,1],I[:,0])
		self.grating.angles = (angleI+np.arctan2(S[:,1],S[:,0]))/2
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
		#import pdb;pdb.set_trace()
		if self.keff == self.convertToKeff(wavelength):
			return
		self.keff = self.convertToKeff(wavelength)
		if self.eval_r1:
			self.r1 = self.grating.points - self.input.center
			self.alpha = np.arctan2(self.r1[:,1],self.r1[:,0])-self.input.angle
			self.alphaI = (np.arctan2(-self.r1[:,1],-self.r1[:,0])-self.grating.angles)
			self.r1 = np.sqrt(np.sum(self.r1**2,1))
			self.eval_r1 = False
		alpha = self.alpha
		x = self.input.x.reshape(self.input.x.shape+(1,))
		E = self.input.E.reshape(self.input.E.shape+(1,))
		k = self.keff
		r1 = self.r1
		dx = self.input.x[1]-self.input.x[0]
		U = np.sum(E*np.exp(-1j*k*np.sin(alpha)*x),0)
		U = U*(1+np.cos(alpha))*np.exp(-1j*k*r1)*dx/(2*np.sqrt(2*np.pi*r1/k))
		self.grating.E = U
			
	def propagateTo(self,target):
		if target.eval_r2:
			target.eval_r2 = False
			O = target.points.reshape(target.points.shape+(1,)).transpose(1,0,2)
			G = self.grating.points.reshape(self.grating.points.shape+(1,)).transpose(1,2,0)
			R = O - G
			target.r2 = np.sqrt(np.sum(R**2,0))
			target.alphaD = self.grating.angles - np.arctan2(R[1],R[0])
			target.sin_diff = (np.sin(self.alphaI)-np.sin(target.alphaD))/2
			target.cos_sum = (np.cos(self.alphaI)+np.cos(target.alphaD))
		E = self.grating.E
		k = self.keff
		w = self.grating.widths
		r2 = target.r2
		sin_dif = target.sin_diff
		cos_sum = target.cos_sum
		pi = np.pi
		alphaD = target.alphaD
		sinc = np.sinc
		exp = np.exp
		sqrt = np.sqrt
		U = E*sinc(k*w*sin_dif/pi)*w*exp(-1j*k*r2)/sqrt(r2)*cos_sum
		target.E = np.sum(U,1)*sqrt(k/8/pi)

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

	def transmissionAngle(self,wavelength,target,angles):
		self.propagateToGrating(wavelength)
		out = []
		for angle in angles:
	#		import pdb;pdb.set_trace()
			self.setApertureOnRowlandAtAngle(target,angle)
			self.propagateTo(target)
			out.append(self.fractionCoupledInto(target))
		return out

	def setInputReturnOverWavelength(self,wavelength):
		lambdaEff = wavelength/self.neff(wavelength)
		self.Input.angle = np.arcsin(lambdaEff*m/(2*a))

	def setApertureOnRowlandAtAngle(self,aperture,angle):
		aperture.setAngle(angle+np.pi)
		self.setApertureCenterOnRowlandCircle(aperture)
		aperture.makePoints()
	
	def findPeakTransmissionAngleAt(self,wavelength,aperture):
		self.propagateToGrating(wavelength)
		N = su.fwhm(np.abs(self.grating.E))
		a = self.grating.pitch
		m = self.order
		deltaLambdaEff = wavelength/(m*N)/self.neff(wavelength)
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
		#import pdb;pdb.set_trace()
		result = optimize.golden(f,brack = (startAngle,endAngle), tol=1e-7,full_output=True)
		return result 
		
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
		plt.figure()
		plt.plot(np.abs(self.grating.E))
		
	def exportSpectrometerDetails(self,filename):
		fd = open(filename,'w')
		self.write(fd)
		fd.close()

def rectModeFunction(x,width):
	return np.ones(x.shape)	

def cosModeFunction(x,width):
	return np.cos(np.pi*x/width)*np.sqrt(2/width)
