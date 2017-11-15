# -*- coding: utf-8 -*-

import copy
from numpy import double, ndarray, searchsorted, delete, s_, mean
from types import ListType
from scipro import SciPro
from acf import ACF
#from constants import *

class Oscillogram(SciPro):
	'''Oscillogram data'''
	def __init__(self, x = None, y = None, ytype = 'lin'):
		SciPro.__init__(self, x, y, ytype = ytype, xtype = 'lin', dtype=double)

class FROGTrace(SciPro):
    '''FROG trace data'''
    def __init__(self, x = None, y = None, ytype = 'lin', xtype = 'wl'):
        SciPro.__init__(self, x, y, ytype = ytype, xtype = 'lin', dtype=double)
        self.xtype = xtype

    def copy(self):
        return FROGTrace( copy.deepcopy(self.x), copy.deepcopy(self.y), ytype=self.ytype, xtype=self.xtype)

    def split( self, *arguments, **keywords):
        #split trace by the wavelength
        retlist = []
        if len(arguments) > 0:
            retval = self.copy()
            retbuf = self.copy()
            for a in arguments:
                if type(a) is ListType or type(a) is ndarray:
                    for suba in a:
                        ind = searchsorted( retbuf.x[1].T[0], suba)
                        if ind != 0:
                            retval.x = retbuf.x[:, :ind]
                            retval.y = retbuf.y[:ind, :]
                            retlist.append( retval.copy())
                            retbuf.x = delete( retbuf.x, s_[:ind], 1)
                            retbuf.y = delete( retbuf.y, s_[:ind], 0)
                        else:
                            retlist.append( None)
                else:
                    ind = searchsorted( retbuf.x[1].T[0], a)
                    if ind != 0:
                        retval.x = retbuf.x[:, :ind]
                        retval.y = retbuf.y[:ind, :]
                        retlist.append( retval.copy())
                        retbuf.x = delete( retbuf.x, s_[:ind], 1)
                        retbuf.y = delete( retbuf.y, s_[:ind], 0)
                    else:
                        retlist.append( None)
            if retbuf.x.size > 0:
                retlist.append( retbuf)
            else:
                    retlist.append( None)
        else:
            retlist.append( self.copy())
        return retlist

    def acf(self):
        '''Reduce trace to autocorrelation function'''
        return ACF( self.x[0][0], self.y.sum(axis=0))
    
    def acfAtWl(self, wl):
        '''Return ACF at the nearest wavelength'''
        ind = self.getWlIndexNearest(wl)
        return self.acfAtIndex(ind)
    
    def acfAtIndex(self, ind):
        '''Return ACF at the selected wavelength index'''
        return ACF( self.x[0][0], self.y[ind])

    def getWlIndexNearest(self, wl):
        indr = searchsorted( self.x[1].T[0], wl)
        indl = searchsorted( self.x[1].T[0], wl, 'left')
        if abs(self.x[1][indr][0] - wl) < abs(self.x[1][indl][0] - wl):
            return indr
        else:
            return indl

    def noiseAvgSub(self):
        '''Subtract noise averaged by X scale'''
        retval = self.copy()
        retval.y -= mean(self.y[0])
        return retval

    def plot(self, *arguments, **keywords):
        '''fuction to plot self spectr\nplot(ptype = 'lin', xl = 'Wavelength, nm', yl = 'Intensity, a.u.')'''
        from matplotlib import cm
        from pylab import colorbar, pcolormesh, xlabel, ylabel
        if not keywords.has_key( 'xl'):
            keywords['xl'] = 'Time, ps'
        if not keywords.has_key( 'yl'):
            keywords['yl'] = 'Wavelength, nm'
        surf = pcolormesh(self.x[0], self.x[1], self.y, cmap=cm.rainbow)
        colorbar(surf)
        xlabel(keywords['xl'])
        ylabel(keywords['yl'])
        #SciPro.plot(self, *arguments, **keywords)

