# -*- coding: utf-8 -*-

from numpy import array, append, flipud, log10, double, arange
from string import atof, atoi
from array import array as ar
from .. oscillogram import Oscillogram

def freadtrc(filename):
    '''This function read data from file of LeCroy oscilloscope'''
    fp = open(filename, 'rb')
    ytype = 'lin'
	
    d = ar('c')
    d.fromfile( fp, 2) #header
    d = ar('c')
    d.fromfile( fp, 9) #all size
    sizeall = int(d.tostring())
    sizetoread = sizeall
    d = ar('c')
    d.fromfile( fp, 16) #header
    sizetoread -= 16

    d = ar('c')
    d.fromfile( fp, 18) #osc
    sizetoread -= 18
    oschz = d

    d = ar('c')
    d.fromfile( fp, 4) #hz
    hzh = d
    sizetoread -= 4

    d = ar('c')
    d.fromfile( fp, 22) #zeroes
    sizetoread -= 22

    d = ar('H')
    d.fromfile( fp, 1) #hz 2002 size
    dsize = d
    sizetoread -= 2

    d = ar('c')
    d.fromfile( fp, 14) #zeroes
    sizetoread -= 14

    d = ar('c')
    d.fromfile( fp, 16) #name
    oscname = d
    sizetoread -= 16

    d = ar('i')
    d.fromfile( fp, 14) #hz 58814
    hzint = d
    sizetoread -= 4*14

    d = ar('f')
    d.fromfile( fp, 2)
    hzf = d
    sizetoread -= 4*2

    d = ar('f')
    d.fromfile( fp, 1)
    dy = d
    sizetoread -= 4

    d = ar('f')
    d.fromfile( fp, 1)
    yoffset = d
    sizetoread -= 4

    d = ar('f')
    d.fromfile( fp, 3)
    hzf2 = d
    sizetoread -= 4*3

    d = ar('f')
    d.fromfile( fp, 1)
    dt = d
    sizetoread -= 4

    d = ar('f')
    d.fromfile( fp, 2)
    hzf3 = d
    sizetoread -= 4*2

    d = ar('i')
    d.fromfile( fp, 1)
    hzi2 = d
    sizetoread -= 4

    d = ar('f')
    d.fromfile( fp, 1)
    hzf4 = d
    sizetoread -= 4

    d = ar('c')
    d.fromfile( fp, 16*5-12)
    sizetoread -= 16*5-12

    d = ar('I')
    d.fromfile( fp, 6) #time stamp?
    sizetoread -= 6*4

    d = ar('I')
    d.fromfile( fp, 7) #time stamp?
    sizetoread -= 7*4

    d = ar('c')
    d.fromfile( fp,2+28) 
    sizetoread -= 2+28

    d = ar('b')
    d.fromfile( fp, sizetoread)

    datay = array( d.tolist())*dy-yoffset
    datax = arange( len(datay))*dt

    return Oscillogram( datax, datay, ytype=ytype)

