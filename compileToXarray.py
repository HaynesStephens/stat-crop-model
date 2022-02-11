import xarray as xr
import pandas as pd

def compilePhase2(model, climate, shift_coord='T', shift_val='0'):
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
'''
models_phase2 = ['LPJmL', 'CARAIB', 'APSIM-UGOE', 'LPJ-GUESS', 'PROMET']
climate_phase2 = ['AgMERRA', 'AgMERRA', 'AgMERRA', 'AgMERRA', 'ERAI']
for model, climate in zip(models_phase2, climate_phase2):
    compilePhase2(model, climate)
'''

# LOOP THROUGH SCENARIOS
'''
T_shifts = [2, 4, 6]
W_shifts = [-50, -30, -10, 10, 30]

for T_val in T_shifts:
    for model, climate in zip(models_phase2, climate_phase2):
        try:
            compilePhase2(model, climate, shift_coord='T', shift_val=T_val)
            print('{0} - T - {1} loaded.'.format(model, T_val))
        except:
            print('{0} - T - {1} DOES NOT EXIST.'.format(model, T_val))
            continue

for W_val in W_shifts:
    for model, climate in zip(models_phase2, climate_phase2):
        try:
            compilePhase2(model, climate, shift_coord='W', shift_val=W_val)
            print('{0} - W - {1} loaded.'.format(model, W_val))
        except:
            print('{0} - W - {1} DOES NOT EXIST.'.format(model, W_val))
            continue
'''

#######################
# THEN GO ON TO PHASE 3
#######################

def compilePhase3(model, climate, scenario):
    # GATHER NAMING INFORMATION FOR FILES
    phase = 'phase3a' if (scenario == 'obsclim') else 'phase3b'
    load_path = '/project2/ggcmi/AgMIP.output/{0}/{1}/{2}/{3}/mai/'.format(model, phase, climate, scenario)   

    scen_time = {'obsclim':'1901_2016', 
                 'historical':'1850_2014',
                 'picontrol':'1850_2100',
                 'ssp126':'2015_2100',
                 'ssp585':'2015_2100'}

    scen_cond = {'obsclim':'2015soc_default',
                 'historical':'2015soc_default',
                 'picontrol':'2015soc_1850co2',
                 'ssp126':'2015soc_2015co2',
                 'ssp585':'2015soc_2015co2'}

    # LOAD DATASET WITH ALL VARIABLES FOR A GIVEN SCENARIO 
    vars = ['yield', 'plantday', 'matyday', 'soilmoist1m']
    vars_time = dict(zip(vars, ['annual', 'annual', 'annual', 'monthly']))

    def writeFilepath(var):
        scen_ext = '' if (scenario == 'obsclim') else 'w5e5_'
        return '{0}_{1}_{2}_{3}_{4}-mai-noirr_global_{5}_{6}.nc'.format(model.lower(), climate, scen_ext+scenario,
                                                                        scen_cond[scenario], var, 
                                                                        vars_time[var], scen_time[scenario])

    file_paths = [load_path + writeFilepath(var) for var in vars]
    for path in file_paths:
        print(path)

    # def loadArr(file_path, var):
    #     ds = xr.open_dataarray(file_path, decode_times=False).rename(var)
    #     start = scen_time[scenario].split('_')[0]
    #     p = ds.time.size
    #     if vars_time[var] == 'annual':
    #         ds['time'] = pd.date_range(start=start, periods=p, freq='A')
    #     elif vars_time[var] == 'monthly':
    #         ds['time'] = pd.date_range(start=start, periods=p, freq='M')
    #     return ds.sel(lat=slice(48.75, 36.25), lon=slice(-103.8, -80.75))
    
    # ds = xr.merge([loadArr(path, var) for path, var in zip(file_paths, vars)])

    # # SAVE DATASET TO SEPARATE DIRECTORY
    # save_dir = '/project2/moyer/ag_data/stat-mod-ds/{0}/'.format(phase)
    # save_name = '{0}_{1}_{2}_{3}_mai-noirr_cornbelt_{4}.nc'.format(model.lower(), climate, scenario,
    #                                                                     scen_cond[scenario], scen_time[scenario])
    # ds.to_netcdf(save_dir+save_name)

models_phase3 = ['LPJmL', 'ACEA', 'LDNDC', 'PROMET']

climate_phase3a = ['gswp3-w5e5']
scenarios_phase3a = ['obsclim']

climate_phase3b = ['gfdl-esm4', 'ipsl-cm6a-lr', 'mpi-esm1-2-hr', 'mri-esm2-0', 'ukesm1-0-ll']
scenarios_phase3b = ['historical', 'picontrol',  'ssp126', 'ssp585']

for climate in climate_phase3a:
    for scenario in scenarios_phase3a:
        out = compilePhase3('LPJmL', climate, scenario)

for climate in climate_phase3b[:1]:
    for scenario in scenarios_phase3b:
        out = compilePhase3('LPJmL', climate, scenario)
