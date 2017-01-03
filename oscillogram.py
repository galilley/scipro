# -*- coding: utf-8 -*-

from numpy import double
from scipro import SciPro
from spectrum import Spectrum

class Oscillogram(SciPro):
	'''Oscillogram data'''
	def __init__(self, x = None, y = None, ytype = 'lin'):
		SciPro.__init__(self, x, y, ytype = ytype, xtype = 'lin', dtype=double)

	def copy(self):
		return Oscillogram( self.x.copy(), self.y.copy(), ytype=self.ytype)

		

