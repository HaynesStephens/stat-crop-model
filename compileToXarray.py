import xarray as xr
import pandas as pd

phase = 'phase2'
models_phase2 = ['LPJmL', 'CARAIB', 'APSIM-UGOE', 'LPJ-GUESS', 'PROMET']
climate_phase2 = ['AgMERRA', 'AgMERRA', 'AgMERRA', 'AgMERRA', 'ERAI']
model, climate = models_phase2[0], climate_phase2[0]
load_path = '/project2/ggcmi/AgMIP.output/{0}/{1}/maize/A0/'.format(model, phase)

vars = ['yield', 'plant-day', 'maty-day', 'trzpah2o', 'runoff']
shift_coord = 'T'
shift_val = 0
file_paths = [load_path + '{0}/{1}_{2}_fullharm_{0}_mai_global_annual_1980_2010_C360_T{3}_W0_N200_A0.nc4'.format(var, model.lower(), climate.lower(), shift_val) for var in vars]
ds = xr.merge([xr.openg_dataset(name, decode_times=False) for name in file_paths])
ds['time'] = pd.date_range(start='12-31-1980', periods=31, freq='A')
ds = ds.sel(lat=slice(48.75, 36.25), lon=slice(-103.8, -80.75), time=slice('1981', '2010'))



