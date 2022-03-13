# -*- coding: utf-8 -*-

from numpy import array, linspace, meshgrid
from array import array as ar
from ..frogtrace import FROGTrace

def fread(filename):
    '''This function read data from file of Mesaphotonics FROGSCAN software'''
    fd = open(filename, 'rb')
    d = ar('f')
    d.fromfile( fd, 4)
    d.byteswap()
    v1 = d.pop(0)
    v2 = d.pop(0)
    v3 = d.pop(0)
    v4 = d.pop(0)
    d.fromfile( fd, int(v1))
    d.byteswap()
    vwl = array(d.tolist())
    d = ar('f')
    d.fromfile( fd, int(v2))
    d.byteswap()
    fd.close()
    vt = 1e-3*linspace(-0.5*v2/v1*v3*v4, 0.5*v2/v1*v3*v4, v2/v1)
    vz = array(d.tolist())
    vzarr = vz.reshape(-1, int(v2/v1))
    print( "file ", filename)
    print( "v1=", v1)
    print( "v2=", v2)
    print( "v3=", v3)
    print( "v4=", v4)
    print( "size1=", len(vzarr))
    print( "size2=", len(vzarr[0]))
    #return vzarr, vx, v1, v2, v3, v4
    return FROGTrace( array(meshgrid(vt,vwl)), vzarr)

