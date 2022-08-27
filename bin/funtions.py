#!/bin/env python
#----------------------------------------------------------------------------------------------
# Nombre: bin.py
# Autor : Victor Saavedra
# Agust 2022
# This script contains the funtions using in the script peak_period_clm.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from scipy.stats import stats
import sys
from scipy import ndimage, signal
from scipy.signal import butter, lfilter

def centroid_valid(hi, max_cent, min_cent):
	mode="image"
	"""
    #This function indicates if the centroid is valid, that is, the image is not dark
	#Input:
	#'max_cent': It is the maximum valid value estimated for the centroid
	#'min_cent': It is the minimum valid value estimated for the centroid

   #Output:
   #'valid': Indicates if the video or image is valid for processing
   #valid = 1; The image is ok to be processed
   #valid = 0; The image is not suitable for processing 
"""
	if mode == "image":
		t=list(range(256))
		ac3=0
		ac4=0
		for i in range (len(hi)):
			ac3=ac3+t[i]*hi[i]
			ac4=ac4+hi[i]
		zr=ac3/ac4
		centroid=zr
		if centroid<= min_cent or centroid>= max_cent:
			valid = 0;
			#print  u"La imagen no es valida por iluminacion"
		else:
			valid = 1
			#print  u"La imagen es valida"
	return valid
#--------------------------------------------------------------------------------------		
def butter_bandpass(Fc1, Fc2, fs, order=10):
    nyq = 0.5 * fs
    low = Fc1 / nyq
    high = Fc2 / nyq
    b, a = butter(order, [low, high], btype='bandpass')
    return b, a
#---------------------------------------------------------------------------------------
def butter_bandpass_filter(data, Fc1, Fc2, fs, order=10):
    b, a = butter_bandpass(Fc1, Fc2, fs, order=order)
    y = lfilter(b, a, data)
    return y
#---------------------------------------------------------------------------------------
def tukeywin(window_length, alpha=0.5):
    '''The Tukey window, also known as the tapered cosine window, can be regarded as a cosine lobe of width \alpha * N / 2
    that is convolved with a rectangle window of width (1 - \alpha / 2). At \alpha = 1 it becomes rectangular, and
    at \alpha = 0 it becomes a Hann window.
 
    We use the same reference as MATLAB to provide the same results in case users compare a MATLAB output to this function
    output
 
    Reference
    ---------
    http://www.mathworks.com/access/helpdesk/help/toolbox/signal/tukeywin.html
 
    '''
    # Special cases
    if alpha <= 0:
        return np.ones(window_length) #rectangular window
    elif alpha >= 1:
        return np.hanning(window_length)
 
    # Normal case
    x = np.linspace(0, 1, window_length)
    w = np.ones(x.shape)
 
    # first condition 0 <= x < alpha/2
    first_condition = x<alpha/2
    w[first_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[first_condition] - alpha/2) ))
 
    # second condition already taken care of
 
    # third condition 1 - alpha / 2 <= x <= 1
    third_condition = x>=(1 - alpha/2)
    w[third_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[third_condition] - 1 + alpha/2))) 
 
    return w
#------------------------------------------------------------------------------------------------------------
def nextpow2(x):
    """returns the smallest power of two that is greater than or equal to the 
    absolute value of x.
  
    This function is useful for optimizing FFT operations, which are 
    most efficient when sequence length is an exact power of two.
  
    :Example:
  
    .. doctest::
  
        >>> x = [255, 256, 257]
        >>> nextpow2(x)
        array([8, 8, 9])
  
    """
    res = np.ceil(np.log2(x))
    return res.astype('int')  #we want integer values only but ceil gives float

#-------------------------------------------------------------------------------------------------------------

