# -*- coding: utf-8 -*-

import numpy as np
from .. field import Field
from .. constants import LIGHT_SPEED

def field_convert_back(f):
    """
    Convert SciPro field into pyofss's one
    return: numpy array, time_window, central_frequency
    """
    return f.y.conjugate(), f.x[-1]-f.x[0], f.central_freq

def field_convert(d,  tw=None, cf=None):
    """
    Convert pyofss field into SciPro's one
    """
    t = np.linspace(-tw/2., tw/2., d.size, False)
    return Field(t, d.conjugate(), cf=cf)

def fread_field(filename, tw=None, cf=None, cwl=None):
    """
    This function read the field saved by pyofss library
    :param string filename: Name of the file
    :param double tw: full Time Window [ps]
    :param double cf: central (carrier) frequency [THz]
    :param double cwl: central (carrier) wavelength (exclusive with cwl) [nm]
    :return: SciPro Field object
    """
    if (cf is not None) and (cwl is not None):
            raise Exception("There should be only one parameter of carrier frequency.")

    if filename.endswith(".npz"):
        d = np.load(filename)['field']
    elif filename.endswith(".npy"):
        d = np.load(filename)
    else:
        try:
            d = np.load(filename + '.npz')['field']
        except:
            d = np.load(filename + '.npy')

    if tw is None:
        tw = d.size

    if cwl is not None:
        cf = (LIGHT_SPEED/cwl)*1e-3

    return field_convert(d, tw, cf)

