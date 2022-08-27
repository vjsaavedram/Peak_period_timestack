#!/bin/env python
# Modified script adapted from the script created by Cesar Cartagena periodopico_timestack.m
# This script calculates the value of the peak period of the wave from timestack images
"""
 ============================================================================================== 
 ---------------------->>>  PEAK PERIOD ESTIMATE FROM TIMESTACK IMAGES  <<<--------------------
 ============================================================================================== 
 +--------------------------------------------------------------------------------------------+ 
 |                                  Autor  : Victor J. Saavedra M.                            | 
 |                                  E-mail : victor.saavedrame@amigo.edu.co                   | 
 |                                  Date   : Diciembre, 2016                                  | 
 |                                  Country,City  : Medellin, Colombia                        |
 |                                  Modified : Agust, 2022                                  |
 |                                 Universidad Catolica Luis Amigó                            | 
 +--------------------------------------------------------------------------------------------+ 
                                                                             
 ============================================================================================== 
"""
#  0.1. Libraries
# -------------------------
import matplotlib.pyplot as plt
import os, sys, time, calendar
sys.path.append("./bin")
import funtions
from scipy import ndimage, signal
from scipy.signal import butter, lfilter
from PIL import Image
import math
import cv2
import numpy as np
import datetime as dt
import dateutil.rrule as rrule
#
#  0.2. Enter parameters for image analysis
max_cent=120
min_cent = 90
FM = 7.5
labels = ['clm']
#  0.3. Load image
img2 = cv2.imread(r'clm_t31_03_2014-02-13-12-00.png')
img = Image.open(r'clm_t31_03_2014-02-13-12-00.png')
#  0.4. Apply grayscale filter to image
#-----------------------------------------------------------------------------------------------
gray = cv2.cvtColor(img2, cv2.COLOR_BGR2GRAY)
hi,bins = np.histogram(gray,255,[0,255]) 
#  0.5  Calculate the histogram of the image to determine whether the image is suitable for processing
#-----------------------------------------------------------------------------------------------
valid1 = funtions.centroid_valid (hi,max_cent, min_cent)
if valid1== 0:
    print  (u"The image is invalid due to illumination")
else:
    size = 4500, 651
    print  (u"The image is valid")
    rotate = img.rotate(-90, expand=1).resize(size)
    rotate = np.asarray(rotate)
    rotate = cv2.cvtColor(rotate, cv2.COLOR_BGR2GRAY)
    ret2,th2 = cv2.threshold(rotate,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    fil, col = th2.shape[:2]
    peak_period = []
    posinicial=[]
    posfinal = []
    for i in range (0,fil):
        indposinit = np.where(th2[i,:]==0)[0]
        if(len(indposinit) >= col * 0.70):
            if len(posinicial) == 0:
                posinicial.append(i)
        elif len(posfinal)== 0 and len(posinicial)>0 and i > posinicial[0]:
            indposfinal = np.where(th2[i,:] == 255)[0]
            if len(indposfinal) >= col * 0.02:
                posfinal.append(i)
    seriex = []
    nu = 0
    if len(posinicial)>0 and len(posfinal)>0:
        for k in range (posinicial[0], posfinal[0]+1):
            ind = np.where(th2[k,:] == 0)[0]
            if len(ind) >= col-1:
                seriex.append(rotate[k,:])
            nu += 1 
    if len(seriex)==0 or len(seriex)<10:
        print (u"The image is invalid due to illumination")
    else:
        seriec= []
        seriex = np.asarray(seriex)
        #-------------------------------------------------------------------
        # Filter pass band
        #-------------------------------------------------------------------
        N=10
        Fc1 = 0.05
        Fc2 = 0.5
        fs=FM
        nyq = 0.5 * fs
        low = float(Fc1)/nyq
        high = float(Fc2) /nyq
        b, a = signal.butter(1, [low,high], 'bandpass')
        seriec = signal.filtfilt(b, a, seriex.transpose())
        seriec = seriec.transpose()
        colc = seriec.shape[1]
        w = funtions.tukeywin(seriec.shape[1],0.05)
        windows = []
        for ii in range (len(seriec)):
            windows.append(np.transpose(w))
        windows = np.asarray(windows)
        seriex2 = seriec*windows[0:seriec.shape[0],0:seriec.shape[1]]
        N = seriec.shape[1]
        Np=2**(funtions.nextpow2(N))
        SP = []
        SP2 = []
        aaa = []
        mediaSP = []
        seriex2 = signal.filtfilt(b, a, seriex2)
        fref=2.0*((FM/2.0)/N)*(np.arange(0,N/2.0))
        for i in range (0, seriex2.shape[0]):
            try:
                dt=0.1
                fo = np.fft.fft(seriex2[i,:],N)
                fr=np.fft.fftfreq(N,dt)
            except:
                print ("%s" % i)
            F1 = fo[1: int(math.floor(len(fo)/2.0))] * np.conj(fo[1: int(math.floor(len(fo)/2.0))])
            F2 = F1/2.0
            F2[0]=F1[0]/N
            F3 = F2/(2.0*((FM/2.0)))
            SP.append(F3);
            if np.isnan(F3[i])==True:
                print('es nan')
        sp = np.transpose(SP)
        for q in range(len(sp)):
            me = np.mean(np.transpose(sp[q]))
            mediaSP.append(me)
        ind2 = np.where(mediaSP==np.max(mediaSP))[0]
        frecuencia_pico = fref[ind2];
        periodo_pico = 1./frecuencia_pico
        print ("The peak period is: %s" %periodo_pico)