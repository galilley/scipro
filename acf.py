# -*- coding: utf-8 -*-

from scipro import SciPro
from numpy import arange, array, append, int32, double, linspace, select, diff
from numpy.fft import *
from pylab import find

class ACF(SciPro):
	'''AutoCorrelation Func'''
	def __init__(self, x = None, y = None):
		SciPro.__init__(self, x, y, ytype = 'lin', xtype = 'lin', dtype=double)

	def split(self):
		'''split series'''
		indmaxsample = find( select( [self.x > (max(self.x)*0.9+min(self.x)*0.1), self.x < (max(self.x)*0.1+min(self.x)*0.9)], [1, 1]) > 0)
		dind = diff( indmaxsample)
		subind = find( dind > max(dind)/2.)
		#subind = insert( subind, 0, 0)
		#subind = append( subind,indmaxsample.size-1)
		indmax = array([], dtype=int)
		for i in xrange( 1, subind.size):
			indmax = append( indmax, indmaxsample[(subind[i]+subind[i-1])/2])
		#find minlen
		samplen = indmax[-1]
		for i in xrange(1,indmax.size):
			if indmax[i]-indmax[i-1] < samplen:
				samplen = indmax[i] - indmax[i-1]
		#extract samples
		sa = []
		for i in xrange(1,indmax.size, 2):
			sa.append( ACF(self.x[indmax[i]-samplen-1:indmax[i]], 
				self.y[indmax[i]-samplen-1:indmax[i]]))
			sa.append( ACF(self.x.take(int32(linspace(indmax[i]+samplen,indmax[i],samplen+1))),
				self.y.take(int32(linspace(indmax[i]+samplen,indmax[i],samplen+1)))))
		return sa

	def splitavg(self):
		'''split series and average them'''
		fl = self.split()
		dx = arange(fl[0].x.size, dtype=double)
		dy = arange(fl[0].x.size, dtype=double)
		dx.fill(0.)
		dy.fill(0.)
		for ind in fl:
			dx+= ind.x
			dy+= ind.y
		dx/=len(fl)
		dy/=len(fl)
		return ACF(dx,dy)

	def fffilter(self, width = 0.99):
		'''FFT filter of ACF'''
		ffbuf = fft(self.y**0.5)
		ffbuf[int32((len(ffbuf)*(1-width))/2):len(ffbuf)-int32((len(ffbuf)*(1-width))/2)-1]*=0.0
		return ACF(self.x, abs(ifft(ffbuf))**2)
	
	def plot(self, *arguments, **keywords):
		if not keywords.has_key( 'xl'):
			keywords['xl'] = 'Time, a.u.'
		if not keywords.has_key( 'yl'):
			keywords['yl'] = 'Intensity, a.u.'
		if not keywords.has_key( 'ptype'):
			keywords['ptype'] = 'lin'
		SciPro.plot( self, *arguments, **keywords)


