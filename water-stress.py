import xarray as xr
import numpy as np
import pandas as pd

model = 'LPJmL'
var = 'transp'
T = 0
basedir = '/project2/ggcmi/AgMIP.output/{0}/phase2/maize/A0/{1}/'.format(model,var)
winf = '{0}_agmerra_fullharm_{1}_mai_global_annual_1980_2010_C360_T{2}_Winf_N200_A0.nc4'.format(model.lower(), var, T)
winf = xr.open_dataset( basedir + winf , decode_times=False)
w0 = '{0}_agmerra_fullharm_{1}_mai_global_annual_1980_2010_C360_T{2}_W0_N200_A0.nc4'.format(model.lower(), var, T)
w0 = xr.open_dataset( basedir + w0 , decode_times=False)
wstress = xr.merge([w0,winf])
wstress['time'] = pd.date_range(start='12-31-1980', periods=31, freq='A')
wstress = wstress.sel(lat=slice(48.75, 36.25), lon=slice(-103.8, -80.75), time=slice('1981', '2010'))
