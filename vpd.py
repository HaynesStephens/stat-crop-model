import numpy as np
import xarray as xr


def VPD_FAO_2010(T_max, T_min, RH):
    """
    RH [%]
    T_min [°C]
    T_max [°C]
    """

    def svp(T):
        return 0.6108 * np.exp((17.27 * T) / (T + 237.3))  # kPa

    e_s = (svp(T_max) + svp(T_min)) / 2
    VPD = e_s * (1 - (RH / 100))  # kPa
    return VPD

def getPhase3Data(model, clim_var):
    phase3_path = '/project2/ggcmi/AgMIP.input/phase3/ISIMIP3/climate_land_only/climate3b'
    hist_path = '{0}/historical/{1}/*{2}*.nc'.format(phase3_path, model, clim_var)
    hist_da = xr.open_mfdataset(hist_path).sel(lat=slice(48.75, 36.25), lon=slice(-103.8, -80.75),
                                               time=slice('1981', '2100'))

    ssp5_path = '{0}/ssp585/{1}/*{2}*.nc'.format(phase3_path, model, clim_var)
    ssp5_da = xr.open_mfdataset(ssp5_path).sel(lat=slice(48.75, 36.25), lon=slice(-103.8, -80.75),
                                               time=slice('1981', '2100'))

    joint_da = xr.concat([hist_da, ssp5_da], dim='time')
    return joint_da

models = ['GFDL-ESM4', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'UKESM1-0-LL']
# for model in models[2:]:
#     print(model)
#     tasmin = getPhase3Data(model, '_tasmin_').tasmin - 273.15
#     print('tasmin loaded.')
#     tasmax = getPhase3Data(model, '_tasmax_').tasmax - 273.15
#     print('tasmax loaded.')
#     hurs   = getPhase3Data(model, '_hurs_').hurs
#     print('hurs loaded.')
#     vpd    = VPD_FAO_2010(tasmax, tasmin, hurs).rename('vpd')
#     print('vpd loaded.')
#     vpd.to_netcdf('/project2/moyer/ag_data/climate_projections/phase3/{0}_r1i1p1f1_w5e5_ssp585_vpd_MDW_daily_1981_2100.nc'.format(model.lower()))
#     print('vpd saved.')
#     print('\n')

save_dir = '/project2/moyer/ag_data/climate_projections/phase3/'
for model in models:
    for cvar in ['tasmin', 'tasmax', 'hurs']:
        da = getPhase3Data(model, '_{0}_'.format(cvar))[cvar]
        da.to_netcdf('{0}/{1}_r1i1p1f1_w5e5_ssp585_{2}_MDW_daily_1981_2100.nc'.format(save_dir, model.lower(), cvar))

