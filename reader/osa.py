# -*- coding: utf-8 -*-

from numpy import array, append, flipud, log10, double
from string import atof, atoi
from ..spectrum import *

def fread(filename):
	'''This function read data from trace'x file of Optical Spectrum Analyzer'''
	fp = open(filename, 'r')
	xtype = 'wl'
	ytype = 'log'

	resln = 0.0
	lsunt = 0
	
	fstr = fp.readline()
	while not fstr.strip() == '\"[TRACE DATA]\"':
		if fstr.find('\"BASEL\"') == 0:
			ytype = 'lin'
		if fstr.find('\"RESLN\"') == 0:
			resln = atof(fstr.split(',')[1])
		if fstr.find('\"LSUNT\"') == 0:
			lsunt = atoi(fstr.split(',')[1])
		if fstr is None:
			break
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

	if lsunt == 0:
		if resln != 0.0:
			if ytype == 'log':
				datay -= 10*log10(resln)
			else:
				datay /= resln
		else:
			return None

	return Spectrum( datax, datay, ytype=ytype, xtype=xtype)

