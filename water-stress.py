import xarray as xr
import numpy as np
import pandas as pd

model = 'LPJmL'
var = 'transp'
T = 0
basedir = '/project2/ggcmi/AgMIP.output/{0}/phase2/maize/A0/{1}/'.format(model,var)
winf = '{0}_agmerra_fullharm_{1}_mai_global_annual_1980_2010_C360_T{2}_Winf_N200_A0.nc4'.format(model.lower(), var, T)
winf = xr.open_dataarray( basedir + winf )
w0 = '{0}_agmerra_fullharm_{1}_mai_global_annual_1980_2010_C360_T{2}_W0_N200_A0.nc4'.format(model.lower(), var, T)
w0 = xr.open_dataarray( basedir + w0 )
