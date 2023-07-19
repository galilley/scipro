from ..frogtrace import FROGTrace
import numpy as np
from PIL import Image, ExifTags
from  json import loads

def fread(filename):
    '''Reads .tiff femtoeasy MSFrog trace'''
    
    # reading exif dictionary
    image = Image.open(filename)
    exif = { ExifTags.TAGS[k]: v for k, v in image.getexif().items() if k in ExifTags.TAGS }

    wl_num = exif['ImageLength']  
    dt_num = exif['ImageWidth']

    model, model_params = exif['Model'].split(';')
    model_params_json = loads(model_params)

    temporal_calibration = model_params_json['temporalCalibration']
    polynom_coefs = model_params_json['spectralCalibration']['polynome']
    n_coefs = len(polynom_coefs)

    lin = np.arange(0.0, wl_num, 1.0)
    monoms = [coef*lin**(n_coefs - 1 - n) for n, coef in enumerate(polynom_coefs)]
    vwl = np.sum(np.array(monoms), axis=0)

    dt_lower = -temporal_calibration*(dt_num-1)/2 * 1e-3 #psec
    dt_upper = temporal_calibration*(dt_num-1)/2 * 1e-3  #psec    
    vt = np.linspace(dt_lower, dt_upper, dt_num)

    amplitude = np.array(Image.open(filename))
    return FROGTrace(np.array(np.meshgrid(vt, vwl)), amplitude)