# -*- coding: utf-8 -*-

from .scipro import SciPro
from numpy import arange, array, append, int32, double, linspace, select, diff,interp,where,zeros,log10
from scipy import integrate
from numpy.fft import *
from .constants import *
#from oscillogram import Oscillogram

class Oscillogram(SciPro):
	'''Oscillogram data'''
	def __init__(self, x = None, y = None, ytype = 'lin'):
		SciPro.__init__(self, x, y, ytype = ytype, xtype = 'lin', dtype=double)

class Spectrum(SciPro):
	'''Spectrum data'''
	def __init__(self, x = None, y = None, ytype = 'lin', xtype = 'wl'):
		SciPro.__init__(self, x, y, ytype = ytype, xtype = 'lin', dtype=double)
		self.xtype = xtype
	
	def copy(self):
		return Spectrum( self.x.copy(), self.y.copy(), ytype=self.ytype, xtype=self.xtype)

	def fromwl(self, wl = None, linI = None, ytype = 'lin' ):
		'''make spectrum object from [wl, nm] as x axis'''
		self.x = wl
		self.y = linI
		self.xtype = 'wl'
		self.ytype = ytype
		return self

	def fromfreq(self, freq = None, linI = None, ytype = 'lin'):
		'''make spectrum object from [freq, GHz] as x axis'''
		self.x = freq
		self.y = linI
		self.xtype = 'freq'
		self.ytype = ytype
		return self

	def tofreq(self):
		'''convert self x axis [wl, nm] to [freq, THz]'''
		if self.xtype is 'wl':
			x = LIGHT_SPEED/self.x*1e-3
			x = x[::-1]
			y = self.y.copy()[::-1]
		else:
			x = self.x
			y = self.y
		return Spectrum(x, y, ytype = self.ytype, xtype = 'freq')

	def towl(self):
		'''convert self x axis [freq, THz] to [wl, nm]'''
		if self.xtype is 'freq':
			x = LIGHT_SPEED/self.x*1e-3
			x = x[::-1]
			y = self.y.copy()[::-1]
		else:
			x = self.x
			y = self.y
		return Spectrum(x, y, ytype = self.ytype, xtype = 'wl')

	def fwhm(self):
		'''return FWHM [nm] parameter of specter'''
		return SciPro.bandwidth(self, -3)

	def cutnoise(self, lev=-20):
		'''cutnoise(self, lev=-20)'''
		if lev > 0.:
			if self.ytype is 'lin':
				ind = where(self.y >= max(self.y)*lev)
			else:
				ind = where(self.y >= max(self.y)+10.*log10(lev))
		elif lev < 0.:
			if self.ytype is 'lin':
				ind = where(self.y >= max(self.y)*(10.**(lev/10.)))
			else:
				ind = where(self.y >= max(self.y)+lev)
		self.x = self.x[ind]
		self.y = self.y[ind]
		return self

	def parabolic(self, p, x):
		'''parabolic(self, p, x)'''
		#'''parabolic(self, x, p = None)'''
		#if p is None: p = self.p
		return p[0]+p[2]*((x-p[1])**2)

	def lnparabolicfit(self):
		'''lnparabolicfit(self):'''
		p0 = [1., self.lpeak(), 50.]
		parabolic = lambda p, x: p[0]+p[2]*((x-p[1])**2)
		func = lambda p, x, y: (parabolic(p, x)-y)
		self.lnI = log(self.y)
		self.p, suc = optimize.leastsq(func, p0, args=(self.x, -self.y), maxfev=10000)
		self.delta0 = self.p[0]
		self.delta1 = self.p[1]
		self.delta2 = self.p[2]
		return self.p

	def omega2mean(self, x0 = None):
		if x0 is None:
			x0 = self.lmean()
		return integrate.trapz(self.tolin().y*(self.x-x0)**2, self.x)/integrate.trapz(self.tolin().y, self.x)

	def fft(self, fakerange = 1.):
		'''return time domain in ps'''
		if self.xtype is 'wl':
			sp = self.tofreq().reverse().equidistant()
		else:
			sp = self
		if self.ytype is 'lin':
			inds =  where( sp.y >= (max(sp.y)/2.))[0]
		else:
			inds =  where( sp.y >= (max(sp.y)-3))[0]
		ind0 = int32(( inds[0]+inds[-1])/2)
		dnum = int32((( fakerange-1)*len(sp.x)))
		ffdata = array( [], dtype = double)
		if self.ytype is 'lin':
			ffdata = append( sp.y[ind0:], zeros( dnum, dtype = double))
			ffdata = append( ffdata, sp.y[:ind0])
		else:
			ffdata = append( 10.**(sp.y[ind0:]/10.), zeros( dnum, dtype = double))
			ffdata = append( ffdata, 10.**(sp.y[:ind0]/10.))
		ffydata = abs( fftshift( ifft( ffdata**0.5)))**2
		fmin = abs(sp.x[1]-sp.x[0])
		ffxdata = linspace(-0.5/fmin, 0.5/fmin, sp.x.size+dnum)
		return Oscillogram(ffxdata, ffydata, ytype='lin')

	def plot(self, *arguments, **keywords):
		'''fuction to plot self spectr\nplot(ptype = 'lin', xl = 'Wavelength, nm', yl = 'Intensity, a.u.')'''
		if 'xl' not in keywords:
			keywords['xl'] = 'Wavelength, nm'
		if 'yl' not in keywords:
			keywords['yl'] = 'Intensity, a.u.'
		if 'ptype' not in keywords:
			keywords['ptype'] = 'lin'
		if self.xtype is 'wl':
			if 'xl' not in keywords:
				keywords['xl'] = 'Wavelength, nm'
			SciPro.plot(self, *arguments, **keywords)
		else:
			if 'xl' not in keywords:
				keywords['xl'] = 'Frequency, THz'
			SciPro.plot(self, *arguments, **keywords)

