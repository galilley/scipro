# -*- coding: utf-8 -*-

from scipy import integrate, optimize
from numpy import *
from types import *

class RFL:
    def __init__(self):
        self.ls = 1310 #RFL wavelength [nm]
        self.lp = 1115 #YDFL wavelength [nm]
        self.alphas = 1.2 #Stokes optical losses [dB/km]
        self.alphap = 2.0 #YDFL optical losses [dB/km]
        self.gR = 1.3 #Raman gain coefficient [(km*W)^(-1)]
        self.Len = 0.370 #Length of fiber, [km]
        self.gam = 3 #Kerr coefficient
        self.K = 1
        self.delta0 = 1.4336
        self.delta1 = 1309.722
        self.delta2 = 51.5108
        self.deltaHT0 = 1.47
        self.deltaHT1 = 1309.69
        self.deltaHT2 = 43.216

        self.alphas = self.alphas/(10*log10(exp(1)))
        self.alphap = self.alphap/(10*log10(exp(1)))

        self.deltaNL=self.deltaNLlin
        #self.deltaNL=self.deltaNLpow

    def __str__(self):
        print 'Stokes wavelength [nm]:\tls = '+str(self.ls)
        print 'Pump wavelength [nm]:\tlp = '+str(self.lp)
        print 'Stokes optical losses [dB/km]:\talphas = '+str(self.alphas*(10*log10(exp(1))))
        print 'YDFL optical losses [dB/km]:\talphap = '+str(self.alphap*(10*log10(exp(1))))
        print 'Raman gain const [(km*W)^(-1)]:\tgR = '+str(self.gR)
        print 'Length of fiber [km]:\tLen = '+str(self.Len)
        print 'kerr coeficient:\tgamma = '+str(self.gam)
        print 'FBGs parameter:\tdelta0 = '+str(self.delta0)
        print 'FBGs parameter:\tdelta1 = '+str(self.delta1)
        print 'FBGs parameter:\tdelta2 = '+str(self.delta2)

    def delta(self, wl):
        return (self.delta0+self.delta2*(wl-self.delta1)**2)

    def deltaHT(self, wl):
        return (self.deltaHT0+self.deltaHT2*(wl-self.deltaHT1)**2)

    def deltaNLlin(self, inI, K = None):
        if K==None: K = self.K
        return self.gam*self.Len*(K*inI+0.1)

    def deltaNLpow(self, inI, K = None):
        if K==None: K = self.K
        return self.gam*self.Len*(K*((inI**3))+0.11)

    def deltaNLconst(self, inI, K = None):
        if K==None: K = self.K
        return K

    def powYDFL(self, inI, K = None):
        if K==None: K = self.K
        return (self.alphap+self.ls/self.lp*self.gR*inI)*(self.delta0+self.deltaNL(inI, K)/3.+2*self.alphas*self.Len)/(2*self.gR*(1-exp(-1*(self.alphap*self.Len+self.ls/self.lp*self.gR*inI*self.Len))));

    def powYDFL2(self, inI, K = None):
        if K==None: K = self.K
        return (self.alphap+self.ls/self.lp*self.gR*inI)*(self.delta0+self.deltaNL(inI, K)/2.+2*self.alphas*self.Len)/(2*self.gR*(1-exp(-2*(self.alphap*self.Len+self.ls/self.lp*self.gR*inI*self.Len))));

    def inpow(self, pYDFL, K = None):
        if K==None: K = self.K
        return optimize.fsolve(lambda inI: self.powYDFL(inI, K)-pYDFL, pYDFL/2.)

    def inpow2(self, pYDFL, K = None):
        if K==None: K = self.K
        return optimize.fsolve(lambda inI: self.powYDFL2(inI, K)-pYDFL, pYDFL/2.)

    def findKin(self, wl, linI, inI = None, p0 = None):
        if inI == None: inI = self.inI
        if p0 == None: p0 = self.K
        tmpwl = wl[where(linI>=max(linI)/2.)]
        self.delta1 = (tmpwl[1]+tmpwl[-1])/2.
        func = lambda p, wl, linI: (self.inSpectr(wl, inI, p)-linI)
        self.K, suc = optimize.leastsq(func, p0, args=(wl, linI), maxfev=10000)
        return self.K

    def findKout(self, wl, linI, inI = None, p0 = None):
        if inI == None: inI = self.inI
        if p0 == None: p0 = self.K
        tmpwl = wl[where(linI>=max(linI)/2.)]
        self.delta1 = (tmpwl[1]+tmpwl[-1])/2.
        func = lambda p, wl, linI: (self.outSpectr(wl, inI, p)-linI)
        self.K, suc =  optimize.leastsq(func, p0, args=(wl, linI), maxfev=10000)
        return self.K

    def findKpow(self, pYDFL, inI = None, p0 = None):
        if inI == None: inI = self.inI
        if p0 == None: p0 = self.K
        return optimize.fsolve(lambda K: self.powYDFL(inI, K)-pYDFL, pYDFL/2.)

    def findKpow2(self, pYDFL, inI = None, p0 = None):
        if inI == None: inI = self.inI
        if p0 == None: p0 = self.K
        return optimize.fsolve(lambda K: self.powYDFL2(inI, K)-pYDFL, pYDFL/2.)

    def findKincomplex(self, pYDFL, wl, linI, p0 = None):
        if p0 == None: p0 = self.K
        tmpwl = wl[where(linI>=max(linI)/2.)]
        self.delta1 = (tmpwl[1]+tmpwl[-1])/2.
        func = lambda p, wl, linI: (self.inSpectr(wl, self.inpow(pYDFL, p), p)-linI)
        self.K, suc = optimize.leastsq(func, p0, args=(wl, linI), maxfev=10000)
        return self.K

    def outSpectr(self, wl, inI, K = None):
        if K==None: K = self.K
        return self.deltaHT(wl)*inI/2.*sqrt(self.delta2/(2*self.deltaNL(inI, K)))/(cosh(pi*(wl-self.delta1)*sqrt(self.delta2/(2*self.deltaNL(inI, K)))))

    def inSpectr(self, wl, inI, K = None, delta1 = None):
        if K==None: K = self.K
        if delta1==None: delta1 = self.delta1
        return inI*sqrt(self.delta2/(2*self.deltaNL(inI, K)))/(cosh(pi*(wl-delta1)*sqrt(self.delta2/(2.*self.deltaNL(inI, K)))))

    def findinSpectr(self, wl, outI):
        return 2*outI/self.delta(wl)

    def outpow(self, pYDFL, K = None):
        if K==None: K = self.K
        wl = linspace(self.delta1-5, self.delta1+5, 10000)
        if type(pYDFL) is type(ndarray) or ListType:
            if type(K) is not (type(ndarray) or ListType):
                tmpk = K
                K = arange(len(pYDFL), dtype=float)
                K.fill(tmpk)
            pout = array([])
            for i in range(len(pYDFL)):
                pout = append(pout, integrate.trapz(0.5*self.delta(wl)*self.inSpectr(wl, self.inpow(pYDFL[i]), K[i]), wl))
        else:
            pcur = pYDFL
            pout = integrate.trapz(0.5*self.delta(wl)*self.inSpectr(wl, self.inpow(pcur), K), wl)
        return pout

    def infwhm(self, pYDFL, K = None):
        if K==None: K = self.K
        wl = linspace(self.delta1-2, self.delta1+2, 10000)
        if type(pYDFL) is type(ndarray) or ListType:
            if type(K) is not (type(ndarray) or ListType):
                tmpk = K
                K = arange(len(pYDFL), dtype=float)
                K.fill(tmpk)
            fwhm = array([])
            for i in range(len(pYDFL)):
                inI = self.inpow(pYDFL[i])
                insp = self.inSpectr(wl, inI, K[i])
                tmpwl = wl[where(insp>=max(insp)/2.)]
                fwhm = append(fwhm, tmpwl[-1]-tmpwl[1])
        else:
            pcur = pYDFL
            inI = self.inpow(pcur)
            insp = self.inSpectr(wl, inI, K)
            tmpwl = wl[where(insp>=max(insp)/2.)]
            fwhm = tmpwl[-1]-tmpwl[1]
        return fwhm

    def outfwhm(self, pYDFL, K = None):
        if K==None: K = self.K
        wl = linspace(self.delta1-2, self.delta1+2, 10000)
        if type(pYDFL) is type(ndarray) or ListType:
            if type(K) is not (type(ndarray) or ListType):
                tmpk = K
                K = arange(len(pYDFL), dtype=float)
                K.fill(tmpk)
            fwhm = array([])
            for i in range(len(pYDFL)):
                inI = self.inpow(pYDFL[i])
                outsp = self.outSpectr(wl, inI, K[i])
                tmpwl = wl[where(outsp>=max(outsp)/2.)]
                fwhm = append(fwhm, tmpwl[-1]-tmpwl[1])
        else:
            pcur = pYDFL
            inI = self.inpow(pcur)
            outsp = self.outSpectr(wl, inI, K)
            tmpwl = wl[where(outsp>=max(outsp)/2.)]
            fwhm = fwhm, tmpwl[-1]-tmpwl[1]
        return fwhm

