from numpy import *
from scipy import optimize
import gdspy
import scientific_utils as su
import ConfigParser
reload(gdspy)


def loadSpectrometer(filename):
	class Holder:
		def __init__(self):
			pass
	s = Holder()
	cfg = ConfigParser.ConfigParser()
	fp = open('spectrometer.txt','r')
	cfg.readfp(fp)
	fp.close()
	radius = cfg.get('Parameters','radius')
	order = cfg.get('Parameters','order')
#	neff  = cfg.get('Parameters','neff')
#	cfg.get('Parameters','aberrationFreeWavelength')
	exec('s.GratingPoints='+cfg.get('Detail','gratingpoints'))
	exec('s.GratingWidths='+cfg.get('Detail','gratingwidths'))
	exec('s.GratingAngles='+cfg.get('Detail','gratingangles'))
	exec('s.InputPoint='	+cfg.get('Detail','InputPoint'		))
	exec('s.InputAngle='	+cfg.get('Detail','InputAngle'		))
	exec('s.InputWidth='	+cfg.get('Detail','InputWidth'		))
	exec('s.OutputPoints='+cfg.get('Detail','OutputPoints'	))
	exec('s.OutputAngles='+cfg.get('Detail','OutputAngles'	))
	exec('s.OutputWidths='+cfg.get('Detail','OutputWidths'	))
	return s

def unitVector(angles):
	return transpose(array([cos(angles),sin(angles)]))
		
def rotate(points,angle,center):
	#for point in points:
	#	print point
	sa = sin(angle)
	ca = cos(angle)
	return [((point[0] - center[0]) * ca - (point[1] - center[1]) * sa + center[0], (point[0] - center[0]) * sa + (point[1] - center[1]) * ca + center[1]) for point in points]

def rotateSpectrometer(spectrometer,angle,center):
	spectrometer.GratingPoints = rotate(spectrometer.GratingPoints,angle,center)
	spectrometer.GratingAngles = 	spectrometer.GratingAngles+ angle
	spectrometer.InputPoint = 		rotate(spectrometer.InputPoint	,angle,center)
	spectrometer.InputAngle =		spectrometer.InputAngle + angle
	spectrometer.OutputPoints= 		rotate(spectrometer.OutputPoints,angle,center)
	spectrometer.OutputAngles = 	spectrometer.OutputAngles + angle

def moveSpectrometer(spectrometer,vector):
	spectrometer.GratingPoints 	= vector + spectrometer.GratingPoints
	spectrometer.InputPoint		= vector + spectrometer.InputPoint
	spectrometer.OutputPoints	= vector + spectrometer.OutputPoints

def makeDiffractionGrating(cell,layer,spectrometer):
	GratingPoints =	spectrometer.GratingPoints
	GratingWidths =	spectrometer.GratingWidths
	GratingAngles =	spectrometer.GratingAngles
	vy = vstack([cos(GratingAngles),sin(GratingAngles)]).transpose()
	vx = vstack([cos(GratingAngles-pi/2),sin(GratingAngles-pi/2)]).transpose()

	nGrooves = len(GratingAngles)
	for i in range(nGrooves-1):
		points = zeros((5,2))
		points[0] = GratingPoints[ i ]-0.5*GratingWidths[i]*vx[ i ] 
		points[1] = GratingPoints[ i ]+0.5*GratingWidths[i]*vx[ i ] 
		points[2] = GratingPoints[i+1]-0.5*GratingWidths[i]*vx[i+1] 
		points[3] = points[2]-array([10e-6,0])
		points[4] = points[0]-array([10e-6,0])
		cell.add(gdspy.Polygon(layer, points))
	points = zeros((4,2))
	points = GratingPoints[i]+0.5*GratingWidths[i]*array([cos(GratingAngles[i]+pi/2),sin(GratingAngles[i]+pi/2)])
	points = GratingPoints[i]-0.5*GratingWidths[i]*array([cos(GratingAngles[i]+pi/2),sin(GratingAngles[i]+pi/2)])
	points = points[1]-array([10e-6,0])

def makeBraggDiffractionGrating(cell,layer,spectrometer,width,period,nPeriods):
	GratingPoints =	spectrometer.GratingPoints
	GratingWidths =	spectrometer.GratingWidths
	GratingAngles =	spectrometer.GratingAngles
	vy = vstack([cos(GratingAngles),sin(GratingAngles)]).transpose()
	vx = vstack([cos(GratingAngles-pi/2),sin(GratingAngles-pi/2)]).transpose()

	nGrooves = len(GratingAngles)
	for i in range(nGrooves):
		for j in range(nPeriods):
			reference = GratingPoints[i] - vy[i]*period*j
			points = zeros((4,2))
			points[0] = reference + vx[i]*GratingWidths[i]/2
			points[1] = reference - vx[i]*GratingWidths[i]/2
			points[2] = reference - vx[i]*GratingWidths[i]/2 - vy[i]*width
			points[3] = reference + vx[i]*GratingWidths[i]/2 - vy[i]*width
			cell.add(gdspy.Polygon(layer, points))


def makeOutputs(cell,layer,spectrometer,Middle):
	TaperLength = 50e-6
	WaveguideWidth = 450e-9
	RidgeWidth = 1.5e-6
	Slope = 5e-6
	RightRadius = 20e-6
	RightMiddleRadius = 20e-6
	LeftMiddleRadius = 20e-6
	LeftRadius = 20e-6
	RightAngle = pi/5
	RightMiddleAngle = pi/3.5
	LeftMiddleAngle = pi/2.5
	LeftAngle = pi/2.8
		
	OutputPoints = spectrometer.OutputPoints
	OutputAngles = spectrometer.OutputAngles - pi
	OutputWidths = spectrometer.OutputWidths
	nOutputs = len(OutputAngles)
	
	s = sign(Middle-arange(nOutputs))	
	vy = vstack([cos(OutputAngles),sin(OutputAngles)]).transpose()
	vx = vstack([cos(OutputAngles-pi/2),sin(OutputAngles-pi/2)]).transpose()
	
	TaperedPoints = OutputPoints + vy*TaperLength
	
	# Side pieces
	for i in range(-1,1):
		points = zeros((4,2))
		points[0] = OutputPoints[i]+s[i]*0.5* OutputWidths[i]*vx[i]
		points[1] = OutputPoints[i]+s[i]*(0.5*OutputWidths[i]+RidgeWidth)*vx[i] 
		points[2] = TaperedPoints[i] + s[i]*(0.5*WaveguideWidth+RidgeWidth)*vx[i]
		points[3] = TaperedPoints[i] + s[i]*0.5 *WaveguideWidth*vx[i] 
		cell.add(gdspy.Polygon(layer,points))
	# Middle pieces	
	for i in range(1,nOutputs):
		points = zeros((4,2))
		points[0] = TaperedPoints[ i ]+0.5*WaveguideWidth*vx[ i ]
		points[1] = OutputPoints[ i ]+0.5*OutputWidths[ i ]*vx[ i ]
		points[2] = OutputPoints[i-1]-0.5*OutputWidths[i-1]*vx[i-1]
		points[3] = TaperedPoints[i-1]-0.5*WaveguideWidth*vx[i-1]
		cell.add(gdspy.Polygon(layer,points))
	# Sides
	angle = [RightAngle+(RightMiddleAngle-RightAngle)*i/(ceil(Middle)-1) for i in range(int(ceil(Middle)))]
	angle = angle+[LeftMiddleAngle+(LeftAngle-LeftMiddleAngle)*i/(nOutputs-ceil(Middle)-1) for i in range(nOutputs-int(ceil(Middle)))]
	Radius = [RightRadius+(RightMiddleRadius-RightRadius)*i/(ceil(Middle)-1) for i in range(int(ceil(Middle)))]
	Radius = Radius+[LeftMiddleRadius+(LeftRadius-LeftMiddleRadius)*i/(nOutputs-ceil(Middle)-1) for i in range(nOutputs-int(ceil(Middle)))]
	Radius[int(floor(Middle))] = 70e-6
	angle[int(floor(Middle))] = pi/2.45
	EndAngles = pi/2 + s*angle
	Disp = Slope*(max(abs(arange(nOutputs)-Middle)) - abs(arange(nOutputs)-Middle))
	paths = []
	# Middle
	for i in range(1,nOutputs):
		p0 = TaperedPoints[i-1]+Disp[i-1]*vy[i-1]
		p1 = TaperedPoints[ i ]+Disp[ i ]*vy[ i ]
		c0 = p0 + s[i-1]*vx[i-1]*Radius[i-1]
		c1 = p1 + s[ i ]*vx[ i ]*Radius[ i ]
		ro0 = Radius[i-1]+s[i-1]*(WaveguideWidth/2+RidgeWidth)
		ro1 = Radius[ i ]-s[ i ]*(WaveguideWidth/2+RidgeWidth)
		ri0 = Radius[i-1]+s[i-1]*WaveguideWidth/2
		ri1 = Radius[ i ]-s[ i ]*WaveguideWidth/2
		points = zeros((4,2))
		points[0] = TaperedPoints[ i ]+0.5*WaveguideWidth*vx[ i ]
		points[1] = TaperedPoints[i-1]-0.5*WaveguideWidth*vx[i-1]
		points[2] = p0-0.5*WaveguideWidth*vx[i-1]
		points[3] = p1+0.5*WaveguideWidth*vx[ i ]
		cell.add(gdspy.Polygon(layer,points))	
		ArcIn0 = c0 + ri0*unitVector(linspace(OutputAngles[i-1]+s[i-1]*pi/2,EndAngles[i-1],80))
		Inter = vstack([c0 + ro0*unitVector(EndAngles[i-1]),c1 + ro1*unitVector(EndAngles[i])])
		ArcIn1 = c1 + ri1*unitVector(linspace(EndAngles[i],OutputAngles[ i ]+s[ i ]*pi/2,80))
		cell.add(gdspy.Polygon(layer,vstack([ArcIn0,Inter,ArcIn1])))
		# creating path
		initPoint = c0 + Radius[i-1]*unitVector(EndAngles[i-1])
		paths.append(gdspy.Path(RidgeWidth, initial_point=initPoint, number_of_paths=2, distance=WaveguideWidth+RidgeWidth))
		paths[i-1].point = initPoint
		paths[i-1].direction = EndAngles[i-1]-pi/2*s[i-1]
	for i in range(-1,1):
		points = zeros((4,2))
		points[0] = TaperedPoints[i]+s[i]*0.5*WaveguideWidth*vx[i]
		points[1] = TaperedPoints[i]+s[i]*(0.5*WaveguideWidth+RidgeWidth)*vx[i]
		points[2] = TaperedPoints[i]+s[i]*(0.5*WaveguideWidth+RidgeWidth)*vx[i]+Disp[i]*vy[i]
		points[3] = TaperedPoints[i]+s[i]*0.5*WaveguideWidth*vx[i]+Disp[i]*vy[i]
		cell.add(gdspy.Polygon(layer,points))	
		center = TaperedPoints[i]+Disp[i]*vy[i] +s[i]*vx[i]*Radius[i]	
		ri = Radius[i]-WaveguideWidth/2
		ro = Radius[i]-(WaveguideWidth/2+RidgeWidth)
		cell.add(gdspy.Round(layer,center,ri,ro,OutputAngles[i]+s[i]*pi/2,EndAngles[i]))
		# creating path
		initPoint = center + Radius[i]*unitVector(EndAngles[i])
		if i == -1:
			paths.append(  gdspy.Path(RidgeWidth, initial_point=initPoint, number_of_paths=2, distance=WaveguideWidth+RidgeWidth))
			paths[i].point = initPoint
			paths[i].direction = EndAngles[i]+pi/2
	return paths

def sprawl(layer,paths,LeftTop,Spacing,BendRadius):
	EndPoints = array([[Spacing*i,0.0] for i in range(len(paths)-1,-1,-1)])+ LeftTop
	for i in range(len(paths)):
		DeltaX = EndPoints[i][0]-paths[i].x
		alpha = (pi/2-paths[i].direction)
		d = sign(DeltaX-sign(DeltaX)*BendRadius*(1-cos(alpha)))
		StartAngle = paths[i].direction+d*pi/2
		if abs(DeltaX)>=BendRadius*(1+cos(alpha)):
			paths[i].arc(layer,BendRadius,StartAngle,pi/2)
			paths[i].direction = pi*(1-d)/2
			paths[i].segment(layer,abs(EndPoints[i][0]-paths[i].x)-BendRadius,paths[i].direction)
			paths[i].arc(layer,BendRadius,-pi/2,-paths[i].direction)
		else:
			theta = arccos(1-(d*DeltaX/BendRadius+(1-cos(alpha)))/2)
			EndAngle = (d+1)*pi/2-d*theta
			paths[i].arc(layer,BendRadius,StartAngle,EndAngle)
			paths[i].arc(layer,BendRadius,EndAngle-d*pi,-pi*(d-1)/2)
		paths[i].segment(layer,abs(EndPoints[i][1]-paths[i].y),'+y')
#	return path
#			
#Spectrometer = loadSpectrometer('spectrometer.txt')
#rotateSpectrometer(Spectrometer,45*pi/180,(0,0))
#cell = gdspy.Cell('paths')
#paths = makeOutputs(Spectrometer)
##for path in paths:
##	path.segment(2, 50e-6,path.direction)
#sprawl(paths,[-0.55e-3,0.4e-3])
#cell.add(paths)
#gdspy.gds_view()	
	
	
