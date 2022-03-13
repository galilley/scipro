# -*- coding: utf-8 -*-

from numpy import array, double

def fread(filename):
	'''This function read data from Irtac .dat files'''
	fp = open(filename, 'r')
	
	fstr = fp.readline()	#read first row
	while not fstr.strip() == 'Series: IAC':
		fstr = fp.readline()
	
	fstr = fp.readline()	#read title
	
	datay = []
	fstr = fp.readline()
	while fstr:
		if fstr.strip() == '': 
			break
		fstr = fstr.strip().split('\t')
		datay.append(fstr[1])
		fstr = fp.readline()
	datay = array(datay).astype(double)
	
	while not fstr.strip() == 'Series: Position':
		fstr = fp.readline()
	
	fstr = fp.readline()	#read title
	
	datax = []
	fstr = fp.readline()
	while fstr:
		if fstr.strip() == '': 
			break
		fstr = fstr.strip().split('\t')
		datax.append(fstr[1])
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
	
	datax = []
	datay = []
	fstr = fp.readline()
	while fstr:
		if fstr.strip() == '': 
			break
		fstr = fstr.strip().split('\t')
		datax.append(fstr[0])
		datay.append(fstr[1])
		fstr = fp.readline()
	fp.close()
	datax = array(datax).astype(double)
	datay = array(datay).astype(double)

	return datax, datay
