#!/usr/bin/env python
import numpy

def angle(vector):
	m=numpy.sqrt(vector[...,0]**2+vector[...,1]**2)
	return numpy.arccos(vector[...,0]/m)*numpy.sign(vector[...,1])

def detectPeak(Y,PeakWidth,NumberOfPeaks=1,X=0):
	x=X
	y=Y
        peaks=numpy.zeros((NumberOfPeaks,2))
	if X==0:
		x=numpy.arange(size(Y))
	for i in range(NumberOfPeaks):
		peakY = y[index]
		peakX = x[index]
		peaks[i,1] = peakX
		peaks[i,0] = peakY
		mask=(x<(peakX-PeakWidth/2))&(x>(peakX+peakWidth/2))
		x=numpy.take(x,mask)
		y=numpy.take(y,mask)

def fwhm(*args):
	if len(args) == 1:
		X = None
		Y = args[0]
	else:
		X = args[0]
		Y = args[1]
	peak = max(Y)
	peakIndex = numpy.argmax(Y)
	leftIndex = numpy.argmin((Y[0:peakIndex]-peak/2)**2)
	rightIndex = numpy.argmin((Y[peakIndex:]-peak/2)**2) 
	if X!=None:
		return X[rightIndex+peakIndex]-X[leftIndex]
	return rightIndex+peakIndex-leftIndex
	
def unitVector(angle):
	return numpy.transpose(numpy.vstack([numpy.cos(angle),numpy.sin(angle)]))
 
def smooth(x,window_len=10,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]
