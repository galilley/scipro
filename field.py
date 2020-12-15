# -*- coding: utf-8 -*-

from numpy import abs,exp,complex64,pi,arctan2,real,imag,linspace,zeros,diff,trapz,sqrt
from scipy import optimize, integrate
from scipy.fftpack import fft,ifft,fftshift,fftfreq
from .scipro import SciPro
from .spectrum import Spectrum
from .oscillogram import Oscillogram
import pylab as pl

class Field(SciPro):
    '''
    Field (complex) data, 
    yform could be:
        complex [yr] (default)
        alg [yr+j*yi], 
        exp [yr*exp(j*yi)]
    param string d: Domain, time or freq
    '''
    def __init__(self, x = None, yr = None, yi = None, yform = None, d='time', cf=0.0):
        self.domain = d
        self.central_freq = cf
        if yr is None and yi is None:
            yi = x[2]
            yr = x[1]
            x = x[0]
            if yform is None:
                yform = 'alg'

        if yform is None:
            yform = 'complex'

        if yform == 'complex':
            SciPro.__init__(self, x, yr, ytype = 'lin', xtype = 'lin', dtype=complex64)
        elif yform == 'alg':
            SciPro.__init__(self, x, yr + 1j*yi, ytype = 'lin', xtype = 'lin', dtype=complex64)
        elif yform == 'exp':
            SciPro.__init__(self, x, yr*exp(1j*yr), ytype = 'lin', xtype = 'lin', dtype=complex64)
        else:
            print('unknown yform')
	
    def copy(self):
        return Field( self.x.copy(), self.y.copy(), d=self.domain, cf=self.central_freq)
	
    def phasemerging( self, gap = 4./3):
        retval = self.copy()
        shift = 0;
        for i in range(1, len(self.y)):
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

    def instfreq(self):
        retval = self.copy()
        retval.y = diff(self.phase().phasemerging().y)/(-2*pi*abs(self.x[-1]-self.x[0])/(len(self.x)-1))
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

    def power(self):
        '''return the power value of the data'''
        return integrate.trapz(real(self.y * self.y.conjugate()), self.x)
	
    def normpower(self,pwr=1.):
        '''normalize power of data to pwr'''
        k = pwr/self.power()
        retval = self.copy()
        retval.y = self.y * sqrt(k)
        return retval
    
    def fft(self, asis=False):
        '''
        Fast Fourier transform
        param bool asis: do not perform normalization (True) or keep total energy (False)
        '''
        if self.domain == 'time':
            retval = self.copy()
            retval.domain = 'freq'
            dt = abs(self.x[-1]-self.x[0])/(len(self.x)-1)
            retval.x = fftshift( fftfreq(self.x.size, dt) + self.central_freq)
            retval.y = fftshift( fft( self.y))
            if not asis:
                ten = self.power()
                retval = retval.normpower(ten)
        else:
            raise Exception("fft can not be applied to a frequency domain")
        return retval
    
    def ifft(self, asis=False):
        '''
        inverse Fast Fourier transform
        param bool asis: do not perform normalization (True) or keep total energy (False)
        '''
        if self.domain == 'freq':
            retval = self.copy()
            retval.domain = 'time'
            fmin = abs(self.x[-1]-self.x[0])/(len(self.x)-1)
            retval.x = linspace(-0.5/fmin, 0.5/fmin, self.x.size, False)
            retval.y = ifft( fftshift(self.y))
            if not asis:
                ten = self.power()
                retval = retval.normpower(ten)
        else:
            raise Exception("ifft can not be applied to a time domain")
        return retval

    def abs2(self, p = 2):
        rvy = real(self.y * self.y.conjugate())
        if self.domain == 'time':
            return Oscillogram(self.x, rvy)
        else:
            return Spectrum(self.x, rvy, xtype='freq')

    def plot(self, *arguments, **keywords):
        '''fuction to plot self spectr\nplot(ptype = 'lin', xl = 'Wavelength, nm', yl = 'Intensity, a.u.')'''
        #TODO ptype = log
        #if not keywords.has_key( 'xl'):
        #    keywords['xl'] = 'Wavelength, nm'
        #if not keywords.has_key( 'yl'):
        #    keywords['yl'] = 'Intensity, a.u.'
        if 'pform' in keywords:
            pform = keywords.pop('pform')
        else:
            pform = 'abs'
        if 'pgap' in keywords:
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
        
        if pform == 'abs':
            ax1.plot( self.x, self.abspower().y, *arguments, **keywords)
            ax1.plot( self.x[0], self.abspower().y[0], *arguments, **keywords)
            ax2.plot( self.x[0], self.phase().y[0], *arguments, **keywords)
            ax2.plot( self.x, self.phase().phasemerging(pgap).y, *arguments, **keywords)
            ax1.set_ylabel('Intensity, |A|**2')
            ax2.set_ylabel('Phase, rad')
            if self.domain == 'time':
                ax1.set_xlabel('Time, ps')
            else:
                ax1.set_xlabel('Frequency, THz')
            pl.sca(ax1)
        elif pform == 'real':
            super(Field, self).plot(self.x, real(self.y), *arguments, **keywords)
        elif pform == 'imag':
            super(Field, self).plot(self.x, imag(self.y), *arguments, **keywords)
        else:
            print('Unknown type '+type+', use \"abs\",\"real\" or \"imag\"')

