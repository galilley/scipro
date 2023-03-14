# -*- coding: utf-8 -*-

import numpy as np
from ..frogtrace import FROGTrace
from ..field import Field


def trace_convert_back(trace_s):
    """
    Convert SciPro FROG trace into pypret's one
    return: MeshData object
    """
    from pypret import MeshData
    trace_freq = trace_s.tofreq()
    return MeshData(
            trace_freq.y.T,
            trace_freq.x[0][0]*1e-12,
            trace_freq.x[1].T[0]*2*np.pi*1e12)


def trace_convert(trace_p):
    """
    Convert pypret FROG trace into SciPro's one
    """
    return FROGTrace(
            np.array(
                np.meshgrid(trace_p.axes[0]*1e12,
                            trace_p.axes[1]/(2*np.pi)*1e-12)),
            trace_p.data.T, xtype='freq').towl()


def field_convert_back(f_s):
    """
    Convert SciPro field into pypret's one
    return: Pulse object
    """
    from pypret import FourierTransform, Pulse

    if f_s.domain == 'freq':
        f_s = f_s.ifft()
    ft = FourierTransform(len(f_s.x), dt=(f_s.x[1]-f_s.x[0])*1e-12)
    pulse = Pulse(ft, f_s.central_freq, unit="f")
    pulse.field = f_s.y.conjugate()
    pulse.spectrum = f_s.fft().y.conjugate()
    return pulse


def field_convert(f_p):
    """
    Convert pypret field into SciPro's one
    """
    return Field(f_p.t*1e12, f_p.field.conjugate(), yform='complex', cf=f_p.w0*1e-12/(2*np.pi))
