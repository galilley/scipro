# -*- coding: utf-8 -*-

import numpy as np
from ..frogtrace import FROGTrace


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
