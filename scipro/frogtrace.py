# -*- coding: utf-8 -*-

import copy
from numpy import double, ndarray, searchsorted, delete, s_, mean
from scipy import interpolate
from .scipro import SciPro
from .acf import ACF
from .constants import LIGHT_SPEED


class FROGTrace(SciPro):
    '''FROG trace data'''
    def __init__(self, x=None, y=None, ytype='lin', xtype='wl'):
        SciPro.__init__(self, x, y, ytype=ytype, xtype='lin', dtype=double)
        self.xtype = xtype

    def copy(self):
        return FROGTrace(copy.deepcopy(self.x), copy.deepcopy(self.y), ytype=self.ytype, xtype=self.xtype)

    def tofreq(self):
        '''convert self x axis [wl, nm] to [freq, THz]'''
        if self.xtype == 'wl':
            x = copy.deepcopy(self.x)
            x[1, :] = LIGHT_SPEED/x[1, :]*1e-3
            x = x[:, ::-1]
            y = self.y.copy()[::-1, :]
        else:
            x = self.x
            y = self.y
        return FROGTrace(x, y, ytype=self.ytype, xtype='freq')

    def towl(self):
        '''convert self x axis [freq, THz] to [wl, nm]'''
        if self.xtype == 'freq':
            x = copy.deepcopy(self.x)
            x[1, :] = LIGHT_SPEED/x[1, :]*1e-3
            x = x[:, ::-1]
            y = self.y.copy()[::-1, :]
        else:
            x = self.x
            y = self.y
        return FROGTrace(x, y, ytype=self.ytype, xtype='wl')

    def split(self, *arguments, **keywords):
        # split trace by the wavelength
        retlist = []
        if len(arguments) > 0:
            retval = self.copy()
            retbuf = self.copy()
            for a in arguments:
                if isinstance(a, list) or isinstance(a, ndarray):
                    for suba in a:
                        ind = searchsorted(retbuf.x[1].T[0], suba)
                        if ind != 0:
                            retval.x = retbuf.x[:, :ind]
                            retval.y = retbuf.y[:ind, :]
                            retlist.append(retval.copy())
                            retbuf.x = delete(retbuf.x, s_[:ind], 1)
                            retbuf.y = delete(retbuf.y, s_[:ind], 0)
                        else:
                            retlist.append(None)
                else:
                    ind = searchsorted(retbuf.x[1].T[0], a)
                    if ind != 0:
                        retval.x = retbuf.x[:, :ind]
                        retval.y = retbuf.y[:ind, :]
                        retlist.append(retval.copy())
                        retbuf.x = delete(retbuf.x, s_[:ind], 1)
                        retbuf.y = delete(retbuf.y, s_[:ind], 0)
                    else:
                        retlist.append(None)
            if retbuf.x.size > 0:
                retlist.append(retbuf)
            else:
                retlist.append(None)
        else:
            retlist.append(self.copy())
        return retlist

    def acf(self):
        '''Reduce trace to autocorrelation function'''
        return ACF(self.x[0][0], self.y.sum(axis=0))
    
    def acfAtWl(self, wl):
        '''Return ACF at the nearest wavelength'''
        ind = self.getWlIndexNearest(wl)
        return self.acfAtIndex(ind)
    
    def acfAtIndex(self, ind):
        '''Return ACF at the selected wavelength index'''
        return ACF(self.x[0][0], self.y[ind])

    def getWlIndexNearest(self, wl):
        indr = searchsorted(self.x[1].T[0], wl)
        indl = searchsorted(self.x[1].T[0], wl, 'left')
        if abs(self.x[1][indr][0] - wl) < abs(self.x[1][indl][0] - wl):
            return indr
        else:
            return indl

    def noiseAvgSub(self, ind=6):
        '''Subtract noise averaged by X scale.
            It works well for a whole FROG trace the first spectral pixels of which
            are masked (see OceanOptics dark pixels description and discussion here:
            https://github.com/ap--/python-seabreeze/issues/88).'''
        retval = self.copy()
        retval.y -= mean(self.y[ind, :])
        return retval

    def plot(self, *arguments, **keywords):
        '''fuction to plot self spectr\nplot(ptype = 'lin', xl = 'Wavelength, nm', yl = 'Intensity, a.u.')'''
        from matplotlib import cm
        from pylab import colorbar, pcolormesh, xlabel, ylabel
        if 'xl' not in keywords:
            keywords['xl'] = 'Time, ps'
        if 'yl' not in keywords:
            if self.xtype == 'wl':
                keywords['yl'] = 'Wavelength, nm'
            else:
                keywords['yl'] = 'Frequency, THz'
        surf = pcolormesh(self.x[0], self.x[1], self.y, cmap=cm.rainbow) #cm.nipy_spectral
        colorbar(surf)
        xlabel(keywords['xl'])
        ylabel(keywords['yl'])
        #SciPro.plot(self, *arguments, **keywords)
