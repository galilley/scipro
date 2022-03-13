# -*- coding: utf-8 -*-

from .optics import spectrum
from numpy import sin, mod, linspace, arange, select, log10, pi

crystalLength = 9.45*10**6
crystalPeriod = 12.539*10**3

#properties of the crystal
A = 4.582
B = 0.099169
C = -0.044432
D = -0.02195

def sinc(x):
    if (x == 0):
        return 1.0
    return sin(pi*x)/x

def nri(lam):
    '''returns refractive index from wavelength'''
    return (A + B/((lam*1e-3)**2 + C) \
         + D*((lam*1e-3))**2)**0.5

def dk(lam3, lam1, lam2):
    """returns argument for sinc()
    n3/lam3 - n1/lam1 - n2/lam2 - 1/L"""
    return (nri(lam3)/lam3 - nri(lam1)/lam1 \
            - nri(lam2)/lam2 - 1.0/crystalPeriod)
            #no pi, cause sinc(x) == sin(pi*x)/(pi*x)

#def mfshg(firstHarmonic = SpectrumComb()):
def mfshg(fh = None): #fh := FirstHarmonic

    "gets a first harmonic SpectrumComb, returns the second harmonic SpectrumComb"
#    if (type(SpectrumComb) != SpectrumComb):
#        print'non SpectrumComb variable received by mfshg()'
#        return False
    #checking SpectrumComb validity
    if fh is None:
        return None
    '''
    if (fh.test() == False):
        print'firstHarmonic is non-valid'
        return False'''

    #initialising output spectrum
    sh = spectrum()
    sh.fromfreq(linspace(2.0*fh.freq[0], 2.0*fh.freq[-1], 2.0*len(fh.freq)-1), arange(2.0*len(fh.freq)-1))
    sh.linI.fill(0.0)
    sh.logI.fill(0.0)
    sh.freqtowl()

    for n in range(len(sh.freq)):   #n = 1..shPoints
        p = 0.0 #variable to store power to add to array

        #0 if n is even
        p += (mod(n,2) == 0) and \
            -3.0*(fh.linI[(n)/2]**2)*(sinc(crystalLength*dk(sh.wl[n], \
                        fh.wl[(n)/2], fh.wl[(n)/2]))**2)\
                  or 0.0

        for j in range( (n>=(len(fh.wl)) and\
                         (n - len(fh.wl)+1)) or 0, n/2 + 1):
            p += 4.0 * fh.linI[j] * fh.linI[n-j] \
              * (sinc(crystalLength*dk(sh.wl[n], \
                        fh.wl[j], \
                        fh.wl[n-j])\
                      )**2)
        sh.linI[n] = p
        sh.logI[n] = select([p < 0, p >= 0],[-200, 10*log10(p+1e-8)])
    return sh

