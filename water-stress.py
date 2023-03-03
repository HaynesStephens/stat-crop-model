import xarray as xr
import numpy as np
import pandas as pd

model = 'LPJmL'
var = 'transp'
T = 0
basedir = '/project2/ggcmi/AgMIP.output/{0}/phase2/maize/A0/{1}/'.format(model,var)

def loadArr(model, var, T, W, N='200'):
    da = '{0}_agmerra_fullharm_{1}_mai_global_annual_1980_2010_C360_T{2}_W{3}_N{4}_A0.nc4'.format(model.lower(), var, T, W, N)
    da = xr.open_dataarray( basedir + da , decode_times=False)
    da['time'] = pd.date_range(start='12-31-1980', periods=31, freq='A')
    da = da.sel(lat=slice(48.75, 36.25), lon=slice(-103.8, -80.75), time=slice('1981', '2010'))
    return da

def calcWStress(model, var, T, N='200'):
    w_rf = loadArr(model, var, T, 0, N=N).rename(var+'_rf')
    w_ir = loadArr(model, var, T, 'inf', N=N).rename(var+'_ir')
    w_stress = xr.merge([w_rf, w_ir])
    w_stress = w_stress.where(w_stress[var+'_ir']>0)
    w_stress[var+'_stress'] = 1 - (w_stress[var+'_rf'] / w_stress[var+'_ir'])
    w_stress = w_stress.expand_dims({'Tshift':[T]})
    return w_stress

w_stress = xr.concat([calcWStress(model, var, Ti, N='200') for Ti in np.arange(0,7,2)], dim='Tshift')

