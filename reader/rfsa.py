# -*- coding: utf-8 -*-

from numpy import array, append, double
from ..spectrum import Spectrum

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
			
	datax = []
	datay = []
	fstr = fp.readline()
	while fstr:
		if fstr.strip() == '': 
			break
		fstr = fstr.strip().split(',')
		datax.append(fstr[0])
		datay.append(fstr[1])
		fstr = fp.readline()
	fp.close()
	datax = array(datax).astype(double)
	datay = array(datay).astype(double)
	
	return Spectrum( datax, datay, ytype=ytype, xtype=xtype)

