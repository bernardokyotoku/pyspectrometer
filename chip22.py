#!/usr/bin/env python
from numpy import *
import gdspy
import gds_spectrometer as gdss
reload(gdss)

def MakeSpectrometerCell():
	Spectrometer = gdss.loadSpectrometer('spectrometer.txt',1e+6)
	gdss.rotateSpectrometer(Spectrometer,(45-3.86)*pi/180,(0,0))
	RightBottom = array([0,0])-Spectrometer.GratingPoints[-1]
	gdss.moveSpectrometer(Spectrometer,RightBottom+array([-5,5]))
	GratingCell = gdspy.Cell('DiffractionGrating')
	gdss.makeBraggDiffractionGrating(GratingCell,layer,Spectrometer,width,period,nPeriods)
	OutputCell = gdspy.Cell('Output')
	paths = gdss.makeOutputs(OutputCell,layer,Spectrometer,	Middle = 36.5)
	gdss.sprawl(layer,paths,LeftOutputPosition,OutputSpacing,BendRadius)
	for path in paths:
		path.segment(layer, TaperLength, '+y', final_width = TaperRidgeWidth, final_distance = TaperWidth + TaperRidgeWidth)
		path.segment(layer, PreTaperLength, '+y')
	OutputCell.add(paths)
	SpectrometerCell = gdspy.Cell('Spectrometer')
	InputPath = gdspy.Path(RidgeWidth, initial_point=Spectrometer.InputPoint[0], number_of_paths=2, distance=Spectrometer.InputWidth+RidgeWidth)
	InputPath.segment(layer,InputTaperLength,Spectrometer.InputAngle + pi,final_distance=WaveguideWidth+RidgeWidth)
	InputPath.arc(layer, BendRadius, Spectrometer.InputAngle -3*pi/2,pi/2)
	InputPath.segment(layer,InputPath.x-SpectrometerInputPositionX-BendRadius,'-x')
	InputPath.arc(layer,BendRadius,pi/2,pi)
	InputPath.segment(layer,InputPath.y-SpectrometerInputPositionY,'-y')
	InputCell = gdspy.Cell('Input')
	InputCell.add(InputPath)
	SpectrometerCell.add(InputCell.elements)
	SpectrometerCell.add(OutputCell.elements)
	SpectrometerCell.add(GratingCell.elements)
	return SpectrometerCell
	

def MakeChipSideEdgeMarks(cell,layer,LeftBottom,BaseHeight,thickness,length):
	cell.add(gdspy.Rectangle(layer,(LeftBottom[0],LeftBottom[1]+BaseHeight[1]-length),(LeftBottom[0]+thickness,LeftBottom[1]+BaseHeight[1])))
	cell.add(gdspy.Rectangle(layer,(LeftBottom[0],LeftBottom[1]),(LeftBottom[0]+thickness,LeftBottom[1]+length)))
	cell.add(gdspy.Rectangle(layer,(LeftBottom[0]+BaseHeight[0]-thickness,LeftBottom[1]),(LeftBottom[0]+BaseHeight[0],LeftBottom[1]+length)))
	cell.add(gdspy.Rectangle(layer,(LeftBottom[0]+BaseHeight[0]-thickness,LeftBottom[1]+BaseHeight[1]-length),(LeftBottom[0]+BaseHeight[0],LeftBottom[1]+BaseHeight[1])))
	
def polishMarks(layer,point,direction):
	big = 10
	medium = 10
	small = 5
	interval = 10
	thick= 1
	p=[]
	p.append(gdspy.Rectangle(layer,(point[0]-medium/2,2*interval+point[1]),(point[0]+medium/2,2*interval+point[1]+thick), datatype=0))
	p.append(gdspy.Rectangle(layer,(point[0]-small/2,interval+point[1]),(point[0]+small/2,interval+point[1]+thick), datatype=0))
	p.append(gdspy.Rectangle(layer,(point[0]-big/2,point[1]),(point[0]+big/2,point[1]+thick), datatype=0))
	p.append(gdspy.Rectangle(layer,(point[0]-small/2,-interval+point[1]),(point[0]+small/2,-interval+point[1]+thick), datatype=0))
	p.append(gdspy.Rectangle(layer,(point[0]-medium/2,-2*interval+point[1]),(point[0]+medium/2,-2*interval+point[1]+thick), datatype=0))
	p.append(gdspy.Rectangle(layer,(point[0]-small/2,-3*interval+point[1]),(point[0]+small/2,-3*interval+point[1]+thick), datatype=0))
	p.append(gdspy.Rectangle(layer,(point[0]-medium/2,-4*interval+point[1]),(point[0]+medium/2,-4*interval+point[1]+thick), datatype=0))
	p.append(gdspy.Rectangle(layer,(point[0]-small/2,-5*interval+point[1]),(point[0]+small/2,-5*interval+point[1]+thick), datatype=0))
	if direction == '+y':pass
	elif direction == '-y':
		for i in p:i.rotate(pi,point)
	elif direction == '+x':
		for i in p:i.rotate(-pi/2,point)
	elif direction == '-x':
		for i in p:i.rotate(pi/2,point)
	return p

# Variables definition
width = 165e-3
period = 330e-3
nPeriods = 6
layer = 2
OutputSpacing = 10
BendRadius = 25
RidgeWidth = 1.5
InputTaperLength = 50
WaveguideWidth = 450e-3
SpectrometerInputPositionX = -995
SpectrometerInputPositionY = 1e+3

ChipWidth = 4e+3
ChipLength = 24e+3
StartPoint = array([SpectrometerInputPositionX-3e+3,-ChipWidth/2])
TaperWidth = 120e-3
TaperRidgeWidth = 4
TaperLength = 50
PreTaperLength = 50
InputOutputShift = 3e+3
LengthBeforeBend = 0.5e+3
NumberOfDevices = 10
NumberOfRef = 5
FieldSize = 1e+3
OutputSpacing = 10
LeftOutputPosition = [-985,2.05e+3]

def makeChipFrame():
	LengthBeforeBend = 0.5e+3
	ChipCell = gdspy.Cell('ChipCell')

	for i in range(NumberOfRef):
		TaperStartPoint = StartPoint+array([-10*NumberOfRef+i*10,150])
		path = gdspy.Path(TaperRidgeWidth, TaperStartPoint, number_of_paths = 2, distance = TaperWidth + TaperRidgeWidth)
		path.segment(layer, PreTaperLength, '+y')
		path.segment(layer, TaperLength, '+y', final_width = RidgeWidth, final_distance = WaveguideWidth + RidgeWidth)
		path.segment(layer, LengthBeforeBend+i*10, '+y')
		path.arc(layer, BendRadius, 0, 0.5*pi)
		path.segment(layer, InputOutputShift-2*BendRadius, '-x')
		path.arc(layer, BendRadius, 1.5*pi, pi)
		path.segment(layer, 2.05e+3-path.y, '+y')
		path.segment(layer, TaperLength, '+y', final_width = TaperRidgeWidth, final_distance = TaperWidth + TaperRidgeWidth)
		path.segment(layer, PreTaperLength, '+y')
		ChipCell.add(path)
	ChipCell.add(gdspy.Rectangle(layer,(-5500,-1600),(-5380,-1480)))
	ChipCell.add(gdspy.Rectangle(layer,(6000,0),(6120,120)))
	MakeChipSideEdgeMarks(ChipCell,layer,(-ChipLength/2,-ChipWidth/2),(ChipLength,ChipWidth+1e+3),10,200)
	ChipCell.add(gdspy.Text(layer,'cornell nanophotonics group' , 30,(-ChipLength/2+0.5e+3,-1200)))
	ChipCell.add(gdspy.Text(layer,'bernardo kyotoku' , 30,(-ChipLength/2+0.5e+3,-1240)))
	ChipCell.add(gdspy.Text(layer,'chip 22' , 30,(-ChipLength/2+0.5e+3,-1280)))
	return ChipCell

def vestedSpectrometer():
	cell = gdspy.Cell('VestedCell')
	TaperStartPoint = array([SpectrometerInputPositionX+InputOutputShift,LeftOutputPosition[1]-ChipWidth+TaperLength+PreTaperLength])
	path = gdspy.Path(TaperRidgeWidth, TaperStartPoint, number_of_paths = 2, distance = TaperWidth + TaperRidgeWidth)
	path.segment(layer, PreTaperLength, '+y')
	cell.add(polishMarks(layer,(path.x+125,path.y),'+y'))
	path.segment(layer, TaperLength, '+y', final_width = RidgeWidth, final_distance = WaveguideWidth + RidgeWidth)
	path.segment(layer, LengthBeforeBend, '+y')
	for r in range(3):
		path.arc(layer, BendRadius, 0, 0.5*pi)
		path.segment(layer, FieldSize-2*BendRadius, '-x')
		path.arc(layer, BendRadius, 1.5*pi, pi)
	path.segment(layer, SpectrometerInputPositionY-path.y, '+y')
	cell.add(polishMarks(layer,(-0.995e+3,2.10e+3),'-y'))
	cell.add(path)
#cell.add(gdspy.CellReference(MakeSpectrometerCell(),origin=(0,0)))
	cell.add(MakeSpectrometerCell().elements)
	MakeChipSideEdgeMarks(cell,layer,(-1e+3,-2e+3),(4e+3,5e+3),10,10)
	return cell
VestedCell = vestedSpectrometer()

#for i in range(NumberOfDevices):
#	TaperStartPoint = StartPoint+array([i*1e+3,0])
#	LengthBeforeBend += 10
#	path = gdspy.Path(TaperRidgeWidth, TaperStartPoint, number_of_paths = 2, distance = TaperWidth + TaperRidgeWidth)
#	path.segment(layer, PreTaperLength, '+y')
#	ChipCell.add(polishMarks(layer,(path.x+125,path.y),'+y'))
#	ChipCell.add(gdspy.Text(layer, '%d'%(i) , 30,(path.x+35,path.y+50), angle = -0.5*pi))
#	path.segment(layer, TaperLength, '+y', final_width = RidgeWidth, final_distance = WaveguideWidth + RidgeWidth)
#	path.segment(layer, LengthBeforeBend, '+y')
#	for r in range(3):
#		path.arc(layer, BendRadius, 0, 0.5*pi)
#		path.segment(layer, FieldSize-2*BendRadius, '-x')
#		path.arc(layer, BendRadius, 1.5*pi, pi)
#	path.segment(layer, SpectrometerInputPositionY-path.y, '+y')
#	ChipCell.add(polishMarks(layer,(-6.995e+3+i*1e+3,2.10e+3),'-y'))
#	ChipCell.add(path)


#print 'area=',SpectrometerCell.area()

#array = gdspy.CellArray(SpectrometerCell, 10, 1, (1e+3,0),origin=(-6e+3,0))
#ChipCell.add(array)
gdspy.gds_view([VestedCell])
ChipCell = makeChipFrame()
fp = open('chip22.gds','w')
gdspy.gds_print(fp,cells=[VestedCell], unit=1, precision=1.0e-3)
fp.close()
fp = open('chip22Frame.gds','w')
gdspy.gds_print(fp,cells=[ChipCell], unit=1, precision=1.0e-3)
fp.close()
