#! /usr/bin/env python

######################################################################
# Original script (gaussian2.py) by
# James Miller-Jones  Amsterdam  27/02/06
#
# Script to insert random point sources into a fits image, and store a
# text file with their positions and strengths.  Based on gaussian.py.
# Difference from gaussian.py is that this will run the program a
# user-specified number of times, rather than just once.  It also asks
# for input fluxes in terms of the rms level rather than the maximum
# point in the image, and searches the FITS header for the beam
# parameters, rather than specifying the header line number, which can
# vary.
# 
# Modified by Hanno Spreeuw for numpy and scipy instead of numarray and
# a slightly more extensive use of the convolution routines.
# Not only to make point sources by convolving "delta functions" with 
# Gaussians corresponding to the clean beam to make point sources,
# but also to make correlated noise by convolving uncorrelated noise
# with a dirty beam. This should be pretty close to the noise in actual
# radio maps. There is still a subtle difference in that the noise
# in radio maps comes from (assumed Gaussian) visibility noise transformed
# by an FFT rather than a FT to the image plane.
# Another aspect is the making of maps that test the deconvolution 
# from the clean beam in the tkp source finder. In order to contstruct
# maps for those tests, we do a double convolution. First, we make
# extended sources (with a Gaussian shape) by convolving "delta
# functions" with a Gaussian. Next, we convolve those with the Gaussian
# corresponding to the clean beam.
######################################################################

import pyfits

import numpy as np
from numpy import random as rn
from scipy.signal import convolve,fftconvolve
from scipy import ndimage
import math
import random
import os
import shutil
import glob

def FileExists(filename):
	"""Test whether the file exists.

	Keyword arguments:
	filename -- the name of the file to be checked

	"""

	return len(glob.glob(filename)) > 0

def gauss(x,y):
	"""Compute the value of an elliptical Gaussian.

	Keyword arguments:
	x -- x-coordinate
	y -- y-coordinate

	"""

	# h=(xdim/2)-1
        h=xdim-1
	# k=(ydim/2)
        k=ydim
# Central co-ordinate (accounting for zero-indexing!), and bearing in
# mind that the central pixel is usually 256,257 or 512,513 in AIPS
# (i.e. y one greater than x)

	bpa = posangle

	a = bmaj/(2*math.sqrt(2*math.log(2)))
	b = bmin/(2*math.sqrt(2*math.log(2)))
# Compute ellipse major and minor axes from FWHM of Gaussian
	theta = bpa*math.pi/180
# Convert PA to radians

	distx = x-h
	disty = y-k

	rotx = distx*math.cos(theta)-disty*math.sin(theta)
	roty = distx*math.sin(theta)+disty*math.cos(theta)
	u = ((rotx/a)**2)+((roty/b)**2)

	z = -0.5*u

	return math.e**z
# This computes an elliptical Gaussian convolving function of unit
# amplitude, centred at (h,k), with major axis bmaj, minor axis bmin,
# position angle bpa

def source_gauss(x,y):
	"""Compute the value of an elliptical Gaussian.

	Keyword arguments:
	x -- x-coordinate
	y -- y-coordinate

	"""

	# h=(xdim/2)-1
        h=xdim-1
	# k=(ydim/2)
        k=ydim
# Central co-ordinate (accounting for zero-indexing!), and bearing in
# mind that the central pixel is usually 256,257 or 512,513 in AIPS
# (i.e. y one greater than x)
# Pick a position angle for the source.
	bpa = -0.5*posangle

	a = 2.*bmaj/(2*math.sqrt(2*math.log(2)))
	b = 0.5*bmin/(2*math.sqrt(2*math.log(2)))
# Compute ellipse major and minor axes from FWHM of Gaussian
	theta = bpa*math.pi/180
# Convert PA to radians

	distx = x-h
	disty = y-k

	rotx = distx*math.cos(theta)-disty*math.sin(theta)
	roty = distx*math.sin(theta)+disty*math.cos(theta)
	u = ((rotx/a)**2)+((roty/b)**2)

	z = -0.5*u

	return math.e**z
# This computes an elliptical Gaussian convolving function of unit
# amplitude, centred at (h,k), with major axis 2.*bmaj, minor axis 0.5*bmin,
# position angle -0.5*bpa


def ratohms(radegs):
	"""Convert RA in decimal degrees format to hours, minutes,
	seconds format.

	Keyword arguments:
	radegs -- RA in degrees format

	Return value:
	ra -- list of 3 values, [hours,minutes,seconds]

	"""

	rah = int(radegs/15)
	ram = int((radegs/15-rah)*60)
	ras = (((radegs/15-rah)*60)-ram)*60
	ra = [rah,ram,ras]
	return ra
# Converts a decimal RA to hms format	

def dectodms(decdegs):
	"""Convert Declination in decimal degrees format to hours, minutes,
	seconds format.

	Keyword arguments:
	decdegs -- Dec. in degrees format

	Return value:
	dec -- list of 4 values, [sign (+/-),degrees,minutes,seconds]

	"""
	
	sign = "+"
	if decdegs<0:
		sign = "-"
		decdegs *= -1
	decd = int(decdegs)
	decm = int((decdegs-decd)*60)
	decs = (((decdegs-decd)*60)-decm)*60
	dec = [sign,decd,decm,decs]
	return dec
# Converts a decimal dec to dms format

currentdir = os.getcwd()

if (os.path.exists('random_sources') == 0):
	os.mkdir('random_sources')

os.chdir('random_sources')

inputfile='PURENOISE5024.FITS'

# Now pull out the necessary information from the input image
hdulist = pyfits.open("%s"%currentdir+"/"+"%s"%inputfile)
header=hdulist[0].header
centrera = hdulist[0].header['crval1']
centredec = hdulist[0].header['crval2']
centrera_hms = ratohms(centrera)
centredec_dms = dectodms(centredec)
centrerapix = hdulist[0].header['crpix1']
centredecpix = hdulist[0].header['crpix2']
print("\nImage co-ordinates:\nRA: ",centrera_hms[0],centrera_hms[1],"%5.2f" % centrera_hms[2],"\nDec: ",centredec_dms[0],centredec_dms[1],centredec_dms[2],"%5.2f" % centredec_dms[3])
# Pull out central RA and Dec and the relevant pixel numbers from the header

raincr = hdulist[0].header['cdelt1']
decincr = hdulist[0].header['cdelt2']
raincr_as = raincr*3600
decincr_as = decincr*3600
# Pull out RA and Dec pixel increments from header

f = -1;
while (f < 0):
	cleanbeam = hdulist[0].header[f]
	if (cleanbeam.find('BMAJ') > 0):
		break
	f = f-1
	
cleanbeam = hdulist[0].header[f]
majoraxis = eval(cleanbeam[20:30])
minoraxis = eval(cleanbeam[38:48])
posangle = eval(cleanbeam[53:60])
# Pull out clean beam parameters from header
print("\nRA pixel increment: %6.2e" % raincr_as,"arcsec\nDec pixel increment: %6.2e" % decincr_as,"arcsec")

bmaj = -1*majoraxis/raincr
bmin = minoraxis/decincr
print("\nRestoring beam:\nMajor axis: %7.4f" % bmaj,"pixels\nMinor axis: %7.4f" % bmin,"pixels\nPA: ",posangle,"degs")
# Convert FWHM of major and minor axes to number of pixels

hdulist.close()

# ydim=hdulist[0].header['NAXIS1']
# xdim=hdulist[0].header['NAXIS2']
ydim=2048
xdim=2048

print("\nImage dimensions: ",xdim," x ",ydim,"\n")
source_beam=np.fromfunction(source_gauss, (2*xdim-1,2*ydim-1))
norm_sb=np.sum(source_beam)
print('np.sum(source_beam)=',norm_sb)
clean_beam=np.fromfunction(gauss, (2*xdim-1,2*ydim-1))
norm_cb=np.sum(clean_beam)
print('np.sum(clean_beam)=',norm_cb)
print('product of the two sums=',norm_sb*norm_cb)
print()
# Make a grid of sources
index_arr=np.indices((xdim,ydim))
condition=(np.mod(index_arr[0],32)==0) & (np.mod(index_arr[1],32)==0)
peak=100.
deltas=np.where(condition,peak,0.)
deltas[0,:]=0.
deltas[:,0]=0.
# In this way we will get 49 point sources with a flux of 100 Jy on a regular grid.
sources=fftconvolve(deltas,source_beam,mode='valid')
if os.path.exists('SOURCES.FITS'):
    os.remove('SOURCES.FITS')
pyfits.writeto('SOURCES.FITS',sources,header)
volume_sources=np.sum(sources)
number_sources=np.sum(deltas)/peak
print('This is the total volume of the sources, the number of sources and the volume per source:',volume_sources,number_sources,volume_sources/number_sources)
print()
convolved_sources=fftconvolve(sources,clean_beam,mode='valid')
if os.path.exists('CONVOLVED_SOURCES.FITS'):
    os.remove('CONVOLVED_SOURCES.FITS')
pyfits.writeto('CONVOLVED_SOURCES.FITS',convolved_sources,header)
# Print the maximum values of convolved_sources.
# The peaks of the islands seem to be around 10.*peak, it seems reasonable to use peak as the threshold.
sci_clip = np.where(convolved_sources > peak, 1, 0)
sci_labels, sci_num = ndimage.label(sci_clip)
maximum_values=ndimage.maximum(convolved_sources,sci_labels,np.arange(sci_num+1)[1:])
# print 'The maximum values are:',maximum_values
print('The mean is:',np.mean(maximum_values))
print('The std is:',np.std(maximum_values))
max_max=np.max(maximum_values)
print('The maximum is:',max_max)
min_max=np.min(maximum_values)
print('The minimum is:',min_max)
print('The difference between min_max and max_max is:',max_max-min_max)
print()
# Import the dirty beam for convolution with the uncorrelated pixels.
hdulist_beam = pyfits.open("%s"%currentdir+"/"+"%s"%'DIRTY_BEAM.FITS')
dirty_beam=hdulist_beam[0].data[0,0,:,:]
print('np.sum(dirty_beam)=',np.sum(dirty_beam))
hdulist_beam.close()

no_iter=1
uncorr_data=np.zeros((no_iter,4))
corr_data=np.zeros((no_iter,4))

i=0
while i<no_iter:
    # Make array of random numbers with sigma=1 and mean=0.
    uncorrelated_pixels=rn.randn(256+xdim-1,256+ydim-1)
    flattened_uncorr=np.ravel(uncorrelated_pixels)
    rms_uncorr=np.std(flattened_uncorr)
    mean_uncorr=np.mean(flattened_uncorr)
    uncorr_max=flattened_uncorr.max()
    uncorr_min=flattened_uncorr.min()
    uncorr_data[i,:]=[rms_uncorr,mean_uncorr,uncorr_max,uncorr_min]
    
    correlated_pixels=fftconvolve(uncorrelated_pixels,dirty_beam,mode='valid')
    flattened_corr=np.ravel(correlated_pixels)
    rms_corr=np.std(flattened_corr)
    # print 'The standard deviation of the correlated pixels is',rms_corr
    mean_corr=np.mean(flattened_corr)
    # print 'The mean of the correlated pixels is',mean_corr
    # print 'The significance of this mean against the hypothesis of zero mean, based on 256*256*30/1024 independent  pixels=',mean_corr/(rms_corr*np.sqrt(256*256*30/1024))
    corr_max=flattened_corr.max()
    corr_min=flattened_corr.min()
    # print 'These are the maximum and minimum of the correlated pixels:',corr_max,corr_min
    # print 'The same values in units of the rms noise:',corr_max/rms_corr,corr_min/rms_corr
    corr_data[i,:]=[rms_corr,mean_corr,corr_max,corr_min]

    test_pixels=correlated_pixels+convolved_sources
    if os.path.exists('TEST_DECONV.FITS'):
        os.remove('TEST_DECONV.FITS')
    pyfits.writeto('TEST_DECONV.FITS',test_pixels,header)
    if np.mod(i+1,100)==0:
        print(i+1)
    i+=1

np.save('uncorr_dat',uncorr_data)
np.save('corr_dat',corr_data)

print('This is the significance (in sigma) of the mean of the means of the uncorrelated noise maps:',np.mean(uncorr_data[:,1])/(np.mean(uncorr_data[:,0])/np.sqrt(uncorrelated_pixels.size*no_iter)))

print('This is the significance (in sigma) of the mean of the means of the correlated noise maps:',np.mean(corr_data[:,1])/(np.mean(corr_data[:,0])/np.sqrt(correlated_pixels.size*no_iter*30./1024.)))
if os.path.exists('CORRELATED_NOISE.FITS'):
    os.remove('CORRELATED_NOISE.FITS')
pyfits.writeto('CORRELATED_NOISE.FITS',correlated_pixels,header)


