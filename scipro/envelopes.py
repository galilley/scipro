# -*- coding: utf-8 -*-

from .scipro import SciPro
from .field import Field
from numpy import array, linspace, exp, log, pi, complex, where, select
from numpy.fft import *

dtnum = 2**6

def gaussianField( freq0, width, lev=0.5, xmax=None, num=None):
	'''gaussianField( freq0, width, lev=0.5, xmax=None, num=None)'''
	if num is None and xmax is None:
		x = linspace(-dtnum*width, dtnum*width, 2*dtnum**2)
	elif xmax is None:
		x = linspace(-num*width, num*width, 2*num**2)
	elif num is None:
		x = linspace(-xmax, xmax, 2*(dtnum**2)/xmax)
	else:
		x = linspace(-xmax, xmax, num)
	if lev > 0.:
		y = exp(log(lev)*(2*x/width)**2)
	elif lev < 0.:
		y = exp(log(10.**(lev/10.))*(2*x/width)**2)
	else:
		return None
	y = y*exp(1j*2*pi*freq0*x)
	return Field(x,y)

def gaussianIntensity(width, lev=0.5, xmax=1., num=1024):
	'''gaussianIntensity(width, lev=0.5, xmax=1., num=1024)'''
	if num is None and xmax is None:
		x = linspace(-dtnum*width, dtnum*width, 2*dtnum**2)
	elif xmax is None:
		x = linspace(-num*width, num*width, 2*num**2)
	elif num is None:
		x = linspace(-xmax, xmax, 2*(dtnum**2)/xmax)
	else:
		x = linspace(-xmax, xmax, num)
	if lev > 0.:
		y = exp(log(lev)*(2*x/width)**2)
	elif lev < 0.:
		y = exp(log(10.**(lev/10.))*(2*x/width)**2)
	else:
		return None
	return SciPro(x,y)

def gaussianFieldChirped( freq0, width, beta2, lev=0.5, xmax=None, num=None):
	'''gaussianFieldChirped( freq0, width, beta2, lev=0.5, xmax=None, num=None)'''
	if num is None and xmax is None:
		x = linspace(-dtnum*width, dtnum*width, 2*dtnum**2)
	elif xmax is None:
		x = linspace(-num*width, num*width, 2*num**2)
	elif num is None:
		x = linspace(-xmax, xmax, 2*(dtnum**2)/xmax)
	else:
		x = linspace(-xmax, xmax, num)
	if lev > 0.:
		y = exp(log(lev)*(2*x/width)**2)
	elif lev < 0.:
		y = exp(log(10.**(lev/10.))*(2*x/width)**2)
	else:
		return None
	y = y*exp(1j*2*pi*freq0*x)
	ffy = fft(y)
	freq = fftfreq( len(y), d=abs(x[1]-x[0]))
	ffy = ffy*exp(1j*beta2/2.*(2*pi*(freq-freq0))**2)
	y = ifft(ffy)
	return Field(x,y)

def gaussianIntensityChirped( width, beta2, lev=0.5, xmax=None, num=None):
	'''gaussianIntensityChirped( width, beta2, lev=0.5, xmax=None, num=None)'''
	if num is None and xmax is None:
		x = linspace(-dtnum*width, dtnum*width, 2*dtnum**2)
	elif xmax is None:
		x = linspace(-num*width, num*width, 2*num**2)
	elif num is None:
		x = linspace(-xmax, xmax, 2*(dtnum**2)/xmax)
	else:
		x = linspace(-xmax, xmax, num)
	if lev > 0.:
		y = exp(log(lev)*(2*x/width)**2)
	elif lev < 0.:
		y = exp(log(10.**(lev/10.))*(2*x/width)**2)
	else:
		return None
	ffy = fft(y**0.5)
	freq = fftfreq( len(y), d=abs(x[1]-x[0]))
	#f0 = freq[where(abs(ffy) == max(abs(ffy)))[0][0]]
	#ffy = ffy*exp(1j*2*pi*beta2/2.*(freq-f0)**2)
	ffy = ffy*exp(1j*beta2/2.*(2*pi*freq)**2)
	y = abs(ifft(ffy))**2
	return SciPro(x,y)

def squareIntensity(width, xmax=1., num=1024):
	'''squareIntensity(width, xmax=1., num=1024)'''
	if num is None and xmax is None:
		x = linspace(-dtnum*width, dtnum*width, 2*dtnum**2)
	elif xmax is None:
		x = linspace(-num*width, num*width, 2*num**2)
	elif num is None:
		x = linspace(-xmax, xmax, 2*(dtnum**2)/xmax)
	else:
		x = linspace(-xmax, xmax, num)
	y = select([x < -width/2., x > width/2.], [0.0, 0.0], 1.0)
	return SciPro(x,y)
