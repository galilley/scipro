#!/usr/bin/python
# -*- coding: utf-8 -*-

workpath = '../080924 FBG for japan'
FBGpath = '../080924 FBG for japan/data'
syspath = '../080924 FBG for japan/python'
fnresult = workpath+'/'+'fbg.csv'

isplotvar = True

#
ftype = '.eps'
nstart = 0
nstop = None

#импорт модулей
from pylab import *
from copy import copy
from numpy import array, append, interp
from string import atof
from fpformat import *
import os
import sys
import re

sys.float_output_precision = 4

os.chdir(FBGpath)
sys.path.append(syspath)

from .optics import spectr
import osa

curdir = '.'
ld = os.listdir(FBGpath)

#составление списка файлов

spfbgr = []
spfbgt = []
fnfbgr = []
fnfbgt = []
for fncurrent in ld:
    m = re.search('^([wW](GDF)|(gdf)1060)\w*_(c|C)((\.CSV)$|(\.csv)$)', fncurrent)
    if m: fnfbgr.append(m.string)
    m = re.search('^([wW](GDF)|(gdf)1060)\w*_(d|D)((\.CSV)$|(\.csv)$)', fncurrent)
    if m: fnfbgt.append(m.string)
del m
fnfbgr.sort()
fnfbgt.sort()
fnfbgr = fnfbgr[nstart:nstop]
fnfbgt = fnfbgt[nstart:nstop]

#чтение спектров
for i in range(len(fnfbgr)):
    spfbgr.append(spectrum(osa.fread(fnfbgr[i])))
for i in range(len(fnfbgt)):
    spfbgt.append(spectrum(osa.fread(fnfbgt[i])))


fresult = open(fnresult, 'w')
fresult.write('SN\tdwl-1\tdwl-3\tdwl-10\tcwl\tT[dB]\tR[%]\n')

for i in range(len(fnfbgr)):
    fresult.write(\
    fnfbgr[i].split('_')[1]+'\t'+\
    repr(spfbgr[i].bandwidth(-1).round(decimals=3))+'\t'+\
    repr(spfbgr[i].bandwidth(-3).round(decimals=3))+'\t'+\
    repr(spfbgr[i].bandwidth(-10).round(decimals=3))+'\t'+\
    repr(spfbgr[i].lpeak().round(decimals=3))+'\t'+\
    repr(spfbgt[i].fbgtranspeakloss().round(decimals=1))+'\t'+\
    repr(((1.0-10**(spfbgt[i].fbgtranspeakloss()/10.0))*100).round(decimals=1))+'\n')
fresult.close()

if isplotvar:
    os.chdir(workpath)
    import Gnuplot
    imgdir = './images'
    gp = Gnuplot.Gnuplot()
    gp('set grid xtics ytics lw 3 dashlength 30.0')
    gp('set xtics 1 nomirror')
    gp('set ytics nomirror')
    gp('set mxtics 4')
    gp('set mytics 4')
    gp('set style data points')
    gp('set style function lines')
    
    gp('set border 15 front linetype -1 linewidth 1.000')
    
    gp('set terminal postscript eps enhanced \
        monochrome blacktext \
        linewidth 1.0 butt \
        palfuncparam 2000,0.003 \
        "Helvatica" 24') 
    
    gp('set linecolor black')
    gp('set pointsize 1.0')
    gp('unset key')
    gp('set xrange [ 1098 : 1102 ]')

    
    gp.xlabel('Wavelength [nm]')
    gp.ylabel('Reflected Power [dBm]')
    for i in range(len(fnfbgr)):
        sn = fnfbgr[i].split('_')[1]
        gp.title(sn+' reflect')
        gp('set output '+'\''+imgdir+'/'+'GDF1060_'+sn+'_r'+ftype+'\'')
        gp.plot(Gnuplot.Data(spfbgr[i].wl, spfbgr[i].logI, with_='p pt 7'))

    gp.ylabel('Transmitted Power [dBm]')
    for i in range(len(fnfbgr)):
        sn = fnfbgt[i].split('_')[1]
        gp.title(sn+' transmit')
        gp('set output '+'\''+imgdir+'/'+'GDF1060_'+sn+'_t'+ftype+'\'')
        gp.plot(Gnuplot.Data(spfbgt[i].wl, spfbgt[i].logI, with_='p pt 7'))
    
    gp('set output')
