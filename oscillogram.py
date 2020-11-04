# -*- coding: utf-8 -*-

from numpy import double, where, int32, append, array, zeros
from numpy.fft import fft,fftshift, fftfreq
from .scipro import SciPro
from .spectrum import Spectrum

class Oscillogram(SciPro):
	'''Oscillogram data'''
	def __init__(self, x = None, y = None, ytype = 'lin'):
		SciPro.__init__(self, x, y, ytype = ytype, xtype = 'lin', dtype=double)

	def copy(self):
		return Oscillogram( self.x.copy(), self.y.copy(), ytype=self.ytype)
	def fft(self, fakerange = 1.):
		'''return spectrum domain in THz if the time domain in ps'''
		osc = self
		dnum = int32((( fakerange-1)*len(osc.x)))
		ffdata = array( [], dtype = double)
		if self.ytype is 'lin':
			ffdata = append( osc.y, zeros( dnum, dtype = double))
		else:
			ffdata = append( 10.**(osc.y/10.), zeros( dnum, dtype = double))
		ffydata = abs( fftshift( fft( ffdata)))**2
		dt = abs(osc.x[1]-osc.x[0])
		ffxdata = fftshift(fftfreq(osc.x.size+dnum, dt))
		return Spectrum(ffxdata, ffydata, xtype='freq', ytype='lin')

		

