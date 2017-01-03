# -*- coding: utf-8 -*-

from numpy import array, append, flipud, log10, double
from string import atof
from ..spectrum import *

def fread(filename):
	'''This function read data from trace'x file of RF-Spectrum Analyzer'''
	fp = open(filename, 'r')
	xtype = 'wl'
	ytype = 'log'
	
	fstr = fp.readline()
	while not fstr.strip() == 'DATA':
		if fstr.startswith('Y Axis Units'):
			if fstr.strip().split(',')[1].startswith('dBm'):
				ytype = 'log'
				print 'log'
			else:
				ytype = 'lin'
				print 'lin'
		fstr = fp.readline()
			
	datax = array([], dtype = double)
	datay = array([], dtype = double)
	fstr = fp.readline()
	while fstr:
		if fstr.strip() == '': 
			break
		fstr = fstr.strip().split(',')
		datax = append(datax, atof(fstr[0]))
		datay = append(datay, atof(fstr[1]))
		fstr = fp.readline()
	
	return Spectrum( datax, datay, ytype=ytype, xtype=xtype)

