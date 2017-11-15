# -*- coding: utf-8 -*-

from numpy import array, append, double
from string import atof

def fread(filename):
	'''This function read data from Irtac .dat files'''
	fp = open(filename, 'r')
	
	fstr = fp.readline()	#read first row
	while not fstr.strip() == 'Series: IAC':
		fstr = fp.readline()
	
	fstr = fp.readline()	#read title
	
	datay = array([], dtype=double)
	fstr = fp.readline()
	while fstr:
		if fstr.strip() == '': 
			break
		fstr = fstr.strip().split('\t')
		datay = append(datay, atof(fstr[1]))
		fstr = fp.readline()
	
	while not fstr.strip() == 'Series: Position':
		fstr = fp.readline()
	
	fstr = fp.readline()	#read title
	
	datax = array([], dtype=double)
	fstr = fp.readline()
	while fstr:
		if fstr.strip() == '': 
			break
		fstr = fstr.strip().split('\t')
		datax = append(datax, atof(fstr[1]))
		fstr = fp.readline()	
	
	fp.close()
	return datax, datay

def freadxy(filename):
	'''This function read data from Irtac xy .dat files'''
	fp = open(filename, 'r')
	
	fstr = fp.readline()	#read first row
	while not fstr.strip() == 'Series: IAC':
		fstr = fp.readline()
	
	fstr = fp.readline()	#read title
	
	datax = array([], dtype=double)
	datay = array([], dtype=double)
	fstr = fp.readline()
	while fstr:
		if fstr.strip() == '': 
			break
		fstr = fstr.strip().split('\t')
		datax = append(datax, atof(fstr[0]))
		datay = append(datay, atof(fstr[1]))
		fstr = fp.readline()
	
	fp.close()
	return datax, datay
