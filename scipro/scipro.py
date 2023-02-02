# -*- coding: utf-8 -*-

from pylab import plot, grid, show, xlabel, ylabel, clf
from numpy import alltrue, array, log10, linspace, ndarray, where, append, arange, insert, delete, searchsorted, int32, double, ones, zeros, concatenate, s_, std, arctan2, imag, real, pi, equal, fromfile
from scipy import integrate, optimize, interp
from numpy.fft import fftshift, fft, ifft

# TODO: thread safe plot


class SciPro(object):
    '''this class allow analyze scientific data'''
    def __init__(self, x=None, y=None, ytype='lin', xtype='lin', dtype=double):
        self.ytype = ytype
        self.xtype = xtype
        self.dtype = dtype
        if x is None and y is None:
            self.x = array([], dtype)
            self.y = array([], dtype)
        elif (type(x) is list or type(x) is ndarray or type(x) is tuple) and (len(x) == 2) and (y is None):
            self.x = x[0]
            self.y = x[1]
            self.dtype = type(x[1])
        else:
            try:
                self.x = x.x
                self.y = x.y
                self.xtype = x.xtype
                self.ytype = x.ytype
                self.dtype = type(x.y)
            except:
                self.x = x
                self.y = y
                self.dtype = type(y)
        # reverse X-axis automatically in case of one dimensional array
        if len(self.x.shape) == 1 and len(self.x) > 1 and self.x[0] > self.x[-1]:
            self.x = self.x[::-1]
            self.y = self.y[::-1]

    def __add__(self, var):
        if type(var) is type(self):
            if self.ytype == 'log':
                a = self.convytype('lin')
            else:
                a = self
            if var.ytype == 'log':
                b = var.convytype('lin')
            else:
                b = var
            if alltrue(a.x == b.x):
                y = a.y+b.y
            else:
                y = a.y+interp(a.x, b.x, b.y, left=0., right=0.)
        elif type(var) is int or float:
            if self.ytype == 'log':
                a = self.convytype('lin')
            else:
                a = self
            y = a.y+var
        else:
            return None
        retval = self.copy()
        retval.x = a.x
        retval.y = y
        retval.xtype = a.xtype
        retval.ytype = a.ytype
        return retval

    def __iadd__(self, var):
        if type(var) is type(self):
            if self.ytype == 'log':
                self.setytype('lin')
            if var.ytype == 'log':
                b = var.convytype('lin')
            else:
                b = var
            if alltrue(self.x == b.x):
                self.y += b.y
            else:
                self.y += interp(self.x, b.x, b.y, left=0., right=0.)
        elif type(var) is int or float:
            if self.ytype == 'log':
                self.setytype('lin')
            self.y += var
        else:
            pass
        return self

    def __sub__(self, var):
        if type(var) is type(self):
            if self.ytype == 'log':
                a = self.convytype( 'lin')
            else:
                a = self
            if var.ytype == 'log':
                b = var.convytype( 'lin')
            else:
                b = var
            if alltrue(a.x == b.x):
                y = a.y-b.y
            else:
                y = a.y-interp(a.x, b.x, b.y, left=0., right=0.)
        elif type(var) is int or float:
            if self.ytype == 'log':
                a = self.convytype('lin')
            else:
                a = self
            y = a.y-var
        else:
            return None
        retval = self.copy()
        retval.x = a.x
        retval.y = y
        retval.xtype = a.xtype
        retval.ytype = a.ytype
        return retval

    def __isub__(self, var):
        if type(var) is type(self):
            if self.ytype == 'log':
                self.setytype('lin')
            if var.ytype == 'log':
                b = var.convytype('lin')
            else:
                b = var
            if alltrue(self.x == b.x):
                self.y -= b.y
            else:
                self.y -= interp(self.x, b.x, b.y, left=0., right=0.)
        elif type(var) is int or float:
            if self.ytype == 'log':
                self.setytype( 'lin')
            self.y -= var
        else:
            pass
        return self

    def __neg__(self):
        if self.ytype == 'log':
            a = self.convytype('lin')
        else:
            a = self
        retval = a.copy()
        retval.y *= -1
        return retval

    def __mul__(self, var):
        if type(var) is type(self):
            if self.ytype == 'lin':
                if equal(self.x, var.x).prod():
                    y = self.y*var.y
                else:
                    y = self.y*interp(self.x, var.x, var.y, left=0., right=0.)
            else:
                if equal(self.x, var.x).prod():
                    y = self.y + var.y
                else:
                    y = self.y+interp(self.x, var.x, var.y, left=0., right=0.)
        elif type(var) is int or float:
            if self.ytype == 'lin':
                y = self.y*var
            else:
                y = self.y+10*log10(var)
        else:
            return None
        retval = self.copy()
        retval.y = y
        return retval

    def __imul__(self, var):
        if type(var) is type(self):
            if self.ytype == 'lin':
                self.y *= interp(self.x, var.x, var.y, left=0., right=0.)
            else:
                self.y += interp(self.x, var.x, var.y, left=0., right=0.)
        elif type(var) is int or float:
            if self.ytype == 'lin':
                self.y *= var
            else:
                self.y += 10*log10(var)
        else:
            pass
        return self

    def __div__(self, var):
        if type(var) is type(self):
            if self.ytype == 'lin':
                y = self.y/interp(self.x, var.x, var.y, left=0., right=0.)
            else:
                y = self.y-interp(self.x, var.x, var.y, left=0., right=0.)
        elif type(var) is int or float:
            if self.ytype == 'lin':
                y = self.y/var
            else:
                y = self.y-10*log10(var)
        else:
            return None
        retval = self.copy()
        retval.y = y
        return retval
        
    def __truediv__(self, var):
        return self.__div__(var)

    def __idiv__(self, var):
        if type(var) is type(self):
            if self.ytype == 'lin':
                self.linI /= interp(self.x, var.x, var.y, left=0., right=0.)
            else:
                self.logI -= interp(self.x, var.x, var.y, left=0., right=0.)
        elif type(var) is int or float:
            if self.ytype == 'lin':
                self.y /= var
            else:
                self.y -= 10*log10(var)
        else:
            pass
        return self

    def __getitem__(self, a):  # TODO
        print(a)
        raise AttributeError

    def value(self, x):
        ind = searchsorted(self.x, x)
        return self.y[ind]

    def abspower(self, p=2):
        retval = self.copy()
        retval.y = abs(retval.y)**p
        return retval

    def phasemerging(self, gap=4./3):
        retval = self.copy()
        shift = 0
        for i in range(1, len(self.y)):
            if self.y[i-1] - self.y[i] > gap*pi:
                shift = shift + 2*pi
            elif self.y[i-1] - self.y[i] < -gap*pi:
                shift = shift - 2*pi
            retval.y[i] += shift
        return retval

    def phase(self):
        retval = self.copy()
        retval.y = arctan2(imag(self.y), real(self.y))
        return retval.phasemerging()

    def copy(self):
        return SciPro(self.x.copy(), self.y.copy(), ytype=self.ytype, xtype=self.xtype, dtype=self.dtype)

    def tolog(self):
        if self.ytype == 'log':
            return self
        y = 10*log10(self.y)
        retval = self.copy()
        retval.y = y
        retval.ytype = 'lin'
        return retval

    def tolin(self):
        if self.ytype == 'lin':
            return self
        y = 10.**(self.y/10.)
        retval = self.copy()
        retval.y = y
        retval.ytype = 'lin'
        return retval

    def equidistant(self, num=None, dnum=None):
        '''function make x axis as equidistance by interpolate'''
        retval = self.copy()
        if num is None and dnum is None:
            retval.x = linspace(self.x[0], self.x[-1], len(self.x))
        elif dnum is not None and num is None:
            num = int(abs(self.x[-1]-self.x[0]))/dnum
            retval.x = linspace(self.x[0], self.x[0]+dfreq*(num-1), num)
        elif dnum is None and num is not None:
            retval.x = linspace(self.x[0], self.x[-1], num)
        else:
            dx = (abs(self.x[-1]-self.x[0])-dnum*(num-1))/2.
            retval.x = linspace(self.x[0]+dx, self.x[-1]-dx, num)
        retval.y = interp(retval.x, self.x, self.y)
        return retval

    def bandwidth(self, lev=-3.):
        '''return bandwidth at specified level'''
        if lev > 0.:
            if self.ytype == 'lin':
                cutlev = max(self.y)*lev
            else:
                cutlev = max(self.y)+10.*log10(lev)
        elif lev < 0.:
            if self.ytype == 'lin':
                cutlev = max(self.y)*(10.**(lev/10.))
            else:
                cutlev = max(self.y)+lev
        else:
            cutlev = 0.0
        inds = where(self.y > cutlev)[0]
        s = interp([cutlev], self.y[inds[0]-1:inds[0]+1], self.x[inds[0]-1:inds[0]+1])[0]
        e = interp([cutlev], self.y[inds[-1]:inds[-1]+2], self.x[inds[-1]:inds[-1]+2])[0]
        return e-s

    def bandwidthleft(self, lev=-3.):
        '''return bandwidth at specified level'''
        if lev > 0.:
            if self.ytype == 'lin':
                tmp = self.x[where(self.y >= max(self.y)*lev)]
            else:
                tmp = self.x[where(self.y >= max(self.y)+10.*log10(lev))]
        elif lev < 0.:
            if self.ytype == 'lin':
                tmp = self.x[where(self.y >= max(self.y)*(10.**(lev/10.)))]
            else:
                tmp = self.x[where(self.y >= max(self.y)+lev)]
        else:
            tmp = self.x
        return self.xPeak()-tmp[1]

    def bandwidthright(self, lev=-3.):
        '''return bandwidth at specified level'''
        if lev > 0.:
            if self.ytype == 'lin':
                tmp = self.x[where(self.y >= max(self.y)*lev)]
            else:
                tmp = self.x[where(self.y >= max(self.y)+10.*log10(lev))]
        elif lev < 0.:
            if self.ytype == 'lin':
                tmp = self.x[where(self.y >= max(self.y)*(10.**(lev/10.)))]
            else:
                tmp = self.x[where(self.y >= max(self.y)+lev)]
        else:
            tmp = self.x
        return tmp[-1]-self.xPeak()

    def setxpeak(self, xpeak=0.0):
        '''set peak x in current data'''
        currentxpeak = self.xPeak()
        self.x -= currentxpeak-xpeak
        return self

    def setxzero(self, xzero=0.0):
        '''set peak x in current data'''
        self.x -= xzero
        return self

    def setxzero2peak(self):
        '''set peak x in current data'''
        self.x -= self.xPeak()
        return self

    def setytype(self, ytype):
        if ytype is self.ytype:
            return self
        elif ytype == 'lin':
            self.y = 10.**(self.y/10.)
            self.ytype = 'lin'
        else:
            self.y = 10.*log10(self.y)
            self.ytype = 'log'
        return self

    def convytype(self, ytype):
        if ytype is self.ytype:
            return self.copy()
        elif ytype == 'lin':
            y = 10.**(self.y/10.)
            ytype = 'lin'
        else:
            y = 10.*log10(self.y)
            ytype = 'log'
        retval = self.copy()
        retval.y = y
        retval.ytype = ytype
        return retval

    def cutMin(self, inds=slice(None)):
        retval = self.copy()
        retval.y = self.y - min(self.y[inds])
        return retval

    def xPeak(self):
        '''return x with peak intensity'''
        return self.x[where(self.y == max(self.y))[0][0]]

    def weightedMean(self):
        '''return mean X weighted by Y'''
        return (self.x*self.y).sum()/self.y.sum()

    def xMin(self):
        '''return x with min intensity'''
        return self.x[where(self.linI == min(self.linI))[0][0]]

    def xMean(self):
        '''return x with intensity in center by fwhm'''
        if self.ytype == 'lin':
            tmp = self.wl[where(self.y >= max(self.y)/2.)]
        else:
            tmp = self.wl[where(self.y >= max(self.y)-3.)]
        return (tmp+tmp)/2.0

    def pPeak(self):
        '''return max value'''
        return self.max()

    def max(self):
        '''return max value'''
        return self.y.max()

    def min(self):
        '''return max value'''
        return self.y.min()

    def power(self):
        '''return power of all data'''
        if self.ytype == 'lin':
            return integrate.trapz(self.y, self.x)
        else:
            return integrate.trapz(10.**(self.y/10.), self.x)

    def every(self, skipnum):
        '''прореживание данных, остаётся только каждая skipnum точка'''
        retval = self.copy()
        ds = len(retval.x)
        ds -= ds % skipnum
        retval.x = retval.x[:ds].reshape(-1, skipnum).T[0]
        retval.y = retval.y[:ds].reshape(-1, skipnum).T[0]
        return retval

    def normalize(self):
        '''normalize data to 1'''
        retval = self.copy()
        retval.y = self.y - self.y.min()
        if self.ytype == 'lin':
            retval.y = self.y/self.y.max()
        else:
            retval.y = self.y - self.y.max()
        return retval

    def normpower(self, pwr=1.):
        '''normalize power of data to x'''
        k = pwr/self.power()
        retval = self.copy()
        if self.ytype == 'lin':
            retval.y = self.y * k
        else:
            retval.y = self.y + 10.*log10(k)
        return retval

    def multconst(self, k):
        '''multiply data to constant k'''
        retval = self.copy()
        if self.ytype == 'lin':
            retval.y = self.y * k
        else:
            retval.y = self.y + 10.*log10(k)
        return retval

    def x2Mean(self):
        x0 = self.xMean()
        if self.ytype == 'lin':
            return integrate.trapz(self.y*(self.x-x0)**2, self.x)/integrate.trapz(self.y, self.x)
        else:
            return integrate.trapz(10.**(self.y/10.)*(self.x-x0)**2, self.x)/integrate.trapz(10.**(self.y/10.), self.x)

    def reverse(self):
        retval = self.copy()
        retval.x = retval.x[::-1]
        retval.y = retval.y[::-1]
        return retval

    def fft(self, fakerange=1.):
        '''direct Fourier transform'''
        if self.ytype == 'lin':
            inds = where(abs(self.y) >= (max(abs(self.y))/2.))[0]
        else:
            inds = where(abs(self.y) >= (max(abs(self.y))-3))[0]
        ind0 = int32((inds[0]+inds[-1])/2)
        dnum = int32(((fakerange-1)*len(self.x)))
        ffdata = array([], dtype=double)
        ffdata = append(self.y[ind0:], zeros(dnum, dtype=double))
        ffdata = append(ffdata, self.y[:ind0])
        ffydata = fftshift(fft(ffdata))
        fmin = abs(self.x[1]-self.x[0])
        ffxdata = linspace(-0.5/fmin, 0.5/fmin, self.x.size+dnum)
        return SciPro(ffxdata, ffydata, ytype='lin')

    def ifft(self, fakerange=1.):
        '''inverse Fourier transform'''
        if self.ytype == 'lin':
            inds = where(abs(self.y) >= (max(abs(self.y))/2.))[0]
        else:
            inds = where(abs(self.y) >= (max(abs(self.y))-3))[0]
        ind0 = int32((inds[0]+inds[-1])/2)
        dnum = int32(((fakerange-1)*len(self.x)))
        ffdata = array([], dtype=double)
        ffdata = append(self.y[ind0:], zeros(dnum, dtype=double))
        ffdata = append(ffdata, self.y[:ind0])
        ffydata = fftshift(ifft(ffdata))
        fmin = abs(self.x[1]-self.x[0])
        ffxdata = linspace(-0.5/fmin, 0.5/fmin, self.x.size+dnum)
        return SciPro(ffxdata, ffydata, ytype='lin')

    def autoCorrelationField(self):
        '''int(E(t)*E*(t-tau),dt)'''
        d = arange(self.x.size*2-1, dtype=double)
        for i in range(1, self.x.size+1):
            d[i-1] = d[d.size-i] = integrate.trapz(self.y[-i:]*self.y[:i].conj(), self.x[:i])
        dx = abs(self.x[1]-self.x[0])
        x = linspace(-dx*self.x.size, dx*self.x.size, self.x.size*2-1)
        return SciPro(x, d)

    def autoCorrelationMichelson(self):
        '''int(|E(t)+E(t-tau)|**2,dt)'''
        d = arange(self.x.size*2-1, dtype=double)
        dx = abs(self.x[1]-self.x[0])
        for i in range(1, self.x.size+1):
            valarr = concatenate((self.y[:-i], self.y[-i:]+self.y[:i], self.y[i:]))
            d[i-1] = d[d.size-i] = integrate.trapz( valarr*valarr.conj(), dx=dx)
        x = linspace(-dx*self.x.size, dx*self.x.size, self.x.size*2-1)
        return SciPro(x, d)

    def autoCorrelationIntensity(self):
        '''int(I(t)*I(t-tau),dt)'''
        d = arange(self.x.size*2-1, dtype=double)
        for i in range(1, self.x.size+1):
            d[i-1] = d[d.size-i] = integrate.trapz(self.y[-i:]*self.y[:i], self.x[:i])
        dx = abs(self.x[1]-self.x[0])
        x = linspace(-dx*self.x.size, dx*self.x.size, self.x.size*2-1)
        return SciPro(x, d)

    def autoCorrelationInterferometric(self):
        '''int(|(E(t)+E(t-tau))**2|**2,dt)'''
        d = arange(self.x.size*2-1, dtype=double)
        dx = abs(self.x[1]-self.x[0])
        for i in range(1, self.x.size+1):
            valarr = concatenate((self.y[:-i]**2, (self.y[-i:]+self.y[:i])**2, self.y[i:]**2))
            d[i-1] = d[d.size-i] = integrate.trapz(valarr*valarr.conj(), dx=dx)
        x = linspace(-dx*self.x.size, dx*self.x.size, self.x.size*2-1)
        return SciPro(x, d)

    def crossCorrelationIntensity(self, val):
        '''int(I1(t)*I2(t-tau),dt)'''
        dx = abs(self.x[-1]-self.x[0])/self.x.size
        tm = max(max(abs(self.x)), max(abs(val.x)))
        d = arange(int(round(tm/dx))*2+1, dtype=double)
        t = linspace(-tm, tm, d.size)
        for i in range(d.size):
            d[i] = integrate.trapz(interp(t, self.x, self.y, 0.0, 0.0)*interp(t + t[-(i+1)], val.x, val.y, 0.0, 0.0), t)
        return SciPro(t, d)

    def fftIntensityFilter(self, width=0.9):
        '''FFT filter'''
        ffbuf = fft(self.y**0.5)
        ffbuf[int32((len(ffbuf)*(1-width))/2):len(ffbuf)-int32((len(ffbuf)*(1-width))/2)-1] *= 0.0
        return SciPro(self.x, abs(ifft(ffbuf))**2)

    def acFilter(self):
        '''ac filter'''
        ffbuf = fft(self.y)
        ffbuf[0] *= 0.0
        return SciPro(self.x, abs(ifft(ffbuf)))

    def smoothing(self, num=50):
        sum = 0.0
        num = int(round(num/2.))*2-1
        if num <= 1:
            return self
        smy = arange(self.x.size, dtype=double)
        for ind in range(self.x.size):
            if ind < int32(num/2):
                sum = sum+self.y[ind]
            elif ind < num:
                sum = sum+self.y[ind]
                smy[ind-int32(num/2)] = sum/(ind+1)
            else:
                sum = sum+self.y[ind]-self.y[ind-num]
                smy[ind-int32(num/2)] = sum/num
        for ind in range(int32(num/2)):
            sum = sum-self.y[self.x.size-num+ind]
            smy[self.x.size-int32(num/2)+ind] = sum/(num-1-ind)
        retval = self.copy()
        retval.y = smy
        return retval

    def movingAvg(self, num):
        return self.smoothing(num)

    def movingStd(self, num):
        num = int(round(num/2.))*2-1
        if num <= 1:
            return self
        mvy = arange(self.x.size, dtype=double)
        for ind in range(self.x.size):
            if ind < int32(num/2):
                mvy[ind] = std(self.y[:ind+int32(num/2)])
            elif ind < self.x.size - int32(num/2):
                mvy[ind] = std(self.y[ind-int32(num/2):ind+int32(num/2)])
            else:
                mvy[ind] = std(self.y[ind-int32(num/2):])
        retval = self.copy()
        retval.y = mvy
        return retval
        
    def split(self, *arguments):
        retlist = []
        if len(arguments) > 0:
            retval = self.copy()
            retbuf = self.copy()
            for a in arguments:
                if type(a) is list or type(a) is ndarray:
                    for suba in a:
                        ind = searchsorted(retbuf.x, suba)
                        if ind != 0:
                            retval.x = retbuf.x[:ind]
                            retval.y = retbuf.y[:ind]
                            retlist.append(retval.copy())
                            retbuf.x = delete(retbuf.x, s_[:ind])
                            retbuf.y = delete(retbuf.y, s_[:ind])
                        else:
                            retlist.append(None)
                else:
                    ind = searchsorted(retbuf.x, a)
                    if ind != 0:
                        retval.x = retbuf.x[:ind]
                        retval.y = retbuf.y[:ind]
                        retlist.append(retval.copy())
                        retbuf.x = delete(retbuf.x, s_[:ind])
                        retbuf.y = delete(retbuf.y, s_[:ind])
                    else:
                        retlist.append(None)
            if retbuf.x.size > 0:
                retlist.append(retbuf)
            else:
                retlist.append(None)
        else:
            retlist.append(self.copy())
        return retlist

    def split_filled(self, *arguments, **keywords):
        if 'lev' in keywords:
            lev = keywords.pop('lev')
        else:
            lev = 1e-15
        retlist = []
        #indstart = 0
        xstart = self.x[0]
        ystart = self.y[0]
        if len(arguments) > 0:
            for a in arguments:
                retval = self.copy()
                #retbuf = self.copy()
                #ind = searchsorted( retbuf.x, a)
                inds = where((self.x >= xstart) & (self.x < a))[0]
                # TODO избавиться от insert и сделать нормально через where
                if inds.size != 0:
                    retval.x = self.x.copy()
                    retval.y = zeros(self.x.size-inds.size, dtype=type(ystart))+lev
                    retval.y = insert(retval.y, ones(inds.size, dtype=int)*inds[0], self.y[inds])
                    retlist.append(retval.copy())
                    #indstart = ind
                    xstart = self.x[inds[-1]]
                else:
                    retlist.append(None)
            if xstart < self.x[-1]:
                retval.x = self.x.copy()
                retval.y = zeros(inds[-1], dtype=type(ystart))+lev
                retval.y = append(retval.y, self.y[inds[-1]:])
                retlist.append(retval)
            else:
                retlist.append(None)
        else:
            retlist.append(self.copy())
        return retlist

    def merging(self, var, p0=1000):
        if type(var) is type(self):
            if len(self.x) < len(var.x):
                num = len(self.x)/4
            else:
                num = len(var.x)/4
            xmean = (self.x[num]-self.x[0])/2.0
            func = lambda p, svar, var: (svar.y[0:num]-interp(svar.x[0:num]+p[0], var.x, var.y, left=0., right=0.))
            self.p, suc = optimize.leastsq(func, [p0], args=(self, var), maxfev=1e8, epsfcn=0.1)
            ind0 = searchsorted(self.x, self.x[0]+xmean)
            ind1 = searchsorted(var.x-self.p[0], self.x[0]+xmean)
            return SciPro(append(var.x[:ind1]-self.p[0], self.x[ind0:]), append(var.y[:ind1], self.y[ind0:]))
        else:
            return None
        
    def plot(self, *arguments, **keywords):
        '''fuction to plot self spectr\nplot(type = 'lin', xl = 'xlabel, a.u.', yl = 'ylabel, a.u.')'''
        if 'xl' in keywords:
            xlabel(keywords.pop('xl'))
        else:
            xlabel('xlabel, a.u.')
        if 'yl' in keywords:
            ylabel(keywords.pop('yl'))
        else:
            ylabel('ylabel, a.u.')
        grid('on')
        if 'ptype' in keywords:
            ptype = keywords.pop('ptype')
        else:
            ptype = 'lin'
        if ptype == 'lin':
            if self.ytype == 'lin':
                plot(self.x, self.y, *arguments, **keywords)
            else:
                plot(self.x, 10.**(self.y/10.), *arguments, **keywords)
        elif ptype == 'log':
            if self.ytype == 'lin':
                plot(self.x, 10.*log10(self.y), *arguments, **keywords)
            else:
                plot(self.x, self.y, *arguments, **keywords)
        else:
            print('Unknown type '+type+', use \"lin\" or \"log\"')
            return False
        return True

    def show(self):
        show()

    def clf(self):
        clf()

    def save(self, filename=None):
        if filename is None:
            print('Can\'t save: undefined filename')
            return False
        d = array([self.x, self.y]).T
        d.tofile(filename+".bin", dtype=double)
        return True

    def open(self, filename=None):
        if filename is None:
            print('Can\'t save: undefined filename')
            return False
        d = fromfile(filename + ".bin").reshape(-1, 2).T
        self.x = d[0]
        self.y = d[1]
        return True

    def reshape(self, shape):
        retval = self.copy()
        retval.x = self.x.reshape(shape)
        retval.y = self.y.reshape(shape)
        return retval
