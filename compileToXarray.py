import xarray as xr

phase = 'phase2'
models_phase2 = ['LPJmL', 'CARAIB', 'APSIM-UGOE', 'LPJ-GUESS', 'PROMET']
climate_phase2 = ['AgMERRA', 'AgMERRA', 'AgMERRA', 'AgMERRA', 'ERAI']
model, climate = models_phase2[0], climate_phase2[0]
load_path = '/project2/ggcmi/AgMIP.output/{0}/{1}/maize/A0/'.format(model, phase)

var = 'yield'
shift_coord = 'T'
shift_range = [0, 2, 4, 6]
file_path = load_path + '{0}_{1}_fullharm_{2}_mai_global_annual_1980_2010_C360_T*_W0_N200_A0.nc4'.format(model.lower(), climate.lower(), var)

ds = xr.openmfdataset(file_path)




