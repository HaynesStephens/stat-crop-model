# IMPORTS
import xarray as xr
import numpy as np
import pandas as pd

def loadArr(model, var, T, W):
    if model == 'CARAIB':
        N = 'NA'
    else:
        N = '200'
    if T==0:
        A=0
    else: 
        A=1
    da = '{0}_agmerra_fullharm_{1}_mai_global_annual_1980_2010_C360_T{2}_W{3}_N{4}_A{5}.nc4'.format(model.lower(), var, T, W, N,A)
    da = xr.open_dataarray( basedir + da , decode_times=False)
    da['time'] = pd.date_range(start='12-31-1980', periods=31, freq='A')
    da = da.sel(lat=slice(48.75, 36.25), lon=slice(-103.8, -80.75), time=slice('1981', '2010'))
    return da

def calcWStress(model, var, T):
    w_rf = loadArr(model, var, T, 0).rename(var+'_rf')
    w_ir = loadArr(model, var, T, 'inf').rename(var+'_ir')
    w_stress = xr.merge([w_rf, w_ir])
    w_stress = w_stress.where(w_stress[var+'_ir']>0)
    w_stress[var+'_stress'] = 1 - (w_stress[var+'_rf'] / w_stress[var+'_ir'])
    w_stress = w_stress.expand_dims({'Tshift':[T]})
    return w_stress

# model = 'LPJmL'
# var = 'transp'
models = ['LPJ-GUESS','pDSSAT', 'GEPIC','PEPIC']# CARAIB, LPJmL, EPIC-TAMU Trust
vars = ['transp']

for model in models:
    # if model in ['LPJmL','CARAIB','LPJ-GUESS','pDSSAT']:
    #     vars = ['etransp', 'transp']
    # else: 
    #     vars = ['etransp']
    for var in vars:
        print(model, var)
        w_stress_list = []
        for Ti in np.arange(0,7,2):
            if Ti==0:
                A=0
            else: 
                A=1
            basedir = '/project2/ggcmi/AgMIP.output/{0}/phase2/maize/A{2}/{1}/'.format(model,var,A)
            try: 
                da = calcWStress(model, var, Ti)
                w_stress_list.append(da)
            except:
                continue
        w_stress = xr.concat(w_stress_list, dim='Tshift')
        w_stress = w_stress.to_dataframe().reset_index()
        w_stress['time'] = w_stress.time.dt.year
        w_stress.to_csv('/project2/moyer/ag_data/wstress/UWA1_{0}_{1}_stress.csv'.format(model,var), index=False)