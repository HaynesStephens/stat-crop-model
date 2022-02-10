import xarray as xr
import pandas as pd

def CompilePhase2(model, climate, shift_coord='T', shift_val='0'):
    phase='phase2'
    load_path = '/project2/ggcmi/AgMIP.output/{0}/{1}/maize/A0/'.format(model, phase)

    # LOAD DATASET WITH ALL VARIABLES FOR A GIVEN SCENARIO 
    vars = ['yield', 'plant-day', 'maty-day', 'trzpah2o']#, 'runoff'] No RUNOFF for PROMET
    scen_vals = {'C':360, 'T':0, 'W':0, 'N':200}
    scen_vals[shift_coord] = shift_val
    if model=='CARAIB': scen_vals['N'] = 'NA'
    shift_str = 'C{0}_T{1}_W{2}_N{3}'.format(scen_vals['C'], scen_vals['T'], scen_vals['W'], scen_vals['N'])

    file_paths = [load_path + '{0}/{1}_{2}_fullharm_{0}_mai_global_annual_1980_2010_{3}_A0.nc4'.format(var, model.lower(), climate.lower(), shift_str) for var in vars]
    ds = xr.merge([xr.open_dataarray(name, decode_times=False).rename(var.split('-')[0]) for name, var in zip(file_paths,vars)])
    ds['time'] = pd.date_range(start='12-31-1980', periods=31, freq='A')
    ds = ds.sel(lat=slice(48.75, 36.25), lon=slice(-103.8, -80.75), time=slice('1981', '2010'))
    ds.maty.attrs["units"] = 'days from plant'

    # SAVE DATASET TO SEPARATE DIRECTORY
    save_dir = '/project2/moyer/ag_data/stat-mod-ds/phase2/'
    save_name = '{0}_{1}_cornbelt_1981_2010_{2}_A0.nc4'.format(model.lower(), climate.lower(), shift_str)
    ds.to_netcdf(save_dir+save_name)

# LOOP THROUGH MODELS
models_phase2 = ['LPJmL', 'CARAIB', 'APSIM-UGOE', 'LPJ-GUESS', 'PROMET']
climate_phase2 = ['AgMERRA', 'AgMERRA', 'AgMERRA', 'AgMERRA', 'ERAI']
for model, climate in zip(models_phase2, climate_phase2):
    CompilePhase2(model, climate)

# LOOP THROUGH SCENARIOS
T_shifts = [2, 4, 6]
W_shifts = [-50, -30, -10, 10, 30]

for T_val in T_shifts:
    for model, climate in zip(models_phase2, climate_phase2):
        try:
            CompilePhase2(model, climate, shift_coord='T', shift_val=T_val)
            print('{0} - T - {1} loaded.'.format(model, T_val))
        except:
            print('{0} - T - {1} DOES NOT EXIST.'.format(model, T_val))
            continue

for W_val in W_shifts:
    for model, climate in zip(models_phase2, climate_phase2):
        try:
            CompilePhase2(model, climate, shift_coord='W', shift_val=W_val)
            print('{0} - W - {1} loaded.'.format(model, W_val))
        except:
            print('{0} - W - {1} DOES NOT EXIST.'.format(model, W_val))
            continue


# THEN GO ON TO PHASE 3



