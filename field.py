# -*- coding: utf-8 -*-

from numpy import exp,complex64,pi,arctan2,real,imag,linspace,zeros
from scipy import optimize
from scipro import SciPro
from spectrum import Spectrum
from oscillogram import Oscillogram
from numpy.fft import fft,ifft,fftshift,fftfreq
import pylab as pl

class Field(SciPro):
    '''Field (complex) data, 
    yform could be:
        complex [yr] (default)
        alg [yr+j*yi], 
        exp [yr*exp(j*yi)]'''
    def __init__(self, x = None, yr = None, yi = None, yform = 'complex'):
        if yr is None and yi is None:
            yi = x[2]
            yr = x[1]
            x = x[0]
            yform = 'alg'
        if yform is 'complex':
            SciPro.__init__(self, x, yr, ytype = 'lin', xtype = 'lin', dtype=complex64)
        elif yform is 'alg':
            SciPro.__init__(self, x, yr + 1j*yi, ytype = 'lin', xtype = 'lin', dtype=complex64)
        elif yform is 'exp':
            SciPro.__init__(self, x, yr*exp(1j*yr), ytype = 'lin', xtype = 'lin', dtype=complex64)
        else:
            print('unknown yform')
	
    def copy(self):
        return Field( self.x.copy(), self.y.copy())
	
    def phasemerging( self, gap = 4./3):
        retval = self.copy()
        shift = 0;
        for i in xrange(1, len(self.y)):
            if self.y[i-1] - self.y[i] > gap*pi:
                shift = shift + 2*pi
            elif self.y[i-1] - self.y[i] < -gap*pi:
                shift = shift - 2*pi
            retval.y[i] += shift
        return retval

    def phase( self):
        retval = self.copy()
        retval.y = arctan2(imag(self.y), real(self.y))
        return retval

    def add_phase(self, ph = [0.]):
        retval = self.copy()
        phase = zeros(self.x.size)
        for i in range(len(ph)):
            phase += self.x**i*ph[i]
        retval.y *= exp(1j*phase)
        return retval
    
    def add_chirp(self, chirp_val):
        retval = self.copy()
        retval.y *= exp(1j*retval.x**2*chirp_val)
        return retval
    
    def chirp(self, pgap = 4./3):
        ffpartph = self.phase().phasemerging(gap = pgap)
        funclinfit = lambda p, x: p[0]+p[1]*x+p[2]*x**2
        func = lambda p, x, y: (funclinfit(p, x)-y)
        p0 = [ ffpartph.y[0], (ffpartph.y[-1]-ffpartph.y[0])/(ffpartph.x[-1]-ffpartph.x[0]), 0.]
        p = optimize.leastsq( 
                func, p0, args=(ffpartph.x, ffpartph.y))
        return p[0][2]
    
    def fft(self):
        '''Fast Fourier transform'''
        retval = self.copy()
        dt = abs(self.x[1]-self.x[0])
        retval.x = fftfreq(self.x.size, dt)
        retval.y = fftshift( fft( self.y))
        return retval
    
    def ifft(self):
        '''inverse Fast Fourier transform'''
        retval = self.copy()
        fmin = abs(self.x[1]-self.x[0])
        retval.x = linspace(-0.5/fmin, 0.5/fmin, self.x.size)
        retval.y = ifft( fftshift(self.y))
        return retval

    def plot(self, *arguments, **keywords):
        '''fuction to plot self spectr\nplot(ptype = 'lin', xl = 'Wavelength, nm', yl = 'Intensity, a.u.')'''
        #if not keywords.has_key( 'xl'):
        #    keywords['xl'] = 'Wavelength, nm'
        #if not keywords.has_key( 'yl'):
        #    keywords['yl'] = 'Intensity, a.u.'
        if keywords.has_key( 'pform'):
            pform = keywords.pop('pform')
        else:
            pform = 'abs'
        if keywords.has_key( 'pgap'):
            pgap = keywords.pop('pgap')
        else:
            pgap = 4./3
        
        if len(pl.gcf().axes) > 0:
            ax1 = pl.gcf().axes[0]
        else:
            ax1 = pl.axes()
        if len(pl.gcf().axes) > 1:
            ax2 = pl.gcf().axes[1]
        else:
            ax2 = pl.twinx()
        
        if pform is 'abs':
            ax1.plot( self.x, self.abspower().y, *arguments, **keywords)
            ax1.plot( self.x[0], self.abspower().y[0], *arguments, **keywords)
            ax2.plot( self.x[0], self.phase().y[0], *arguments, **keywords)
            ax2.plot( self.x, self.phase().phasemerging(pgap).y, *arguments, **keywords)
            ax1.set_ylabel('Intensity, |A|**2')
            ax2.set_ylabel('Phase, rad')
            pl.sca(ax1)
        elif pform == 'real':
            pl.plot(self.x, real(self.y), *arguments, **keywords)
        elif pform == 'imag':
            pl.plot(self.x, imag(self.y), *arguments, **keywords)
        else:
            print 'Unknown type '+type+', use \"abs\",\"real\" or \"imag\"'

