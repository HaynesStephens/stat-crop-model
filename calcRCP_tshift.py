import xarray as xr

dat_dir = '/project2/geos39650/ag_data/climate_projections/lpjml'
da   = xr.open_dataarray('{0}/tas_bced_1960_1999_hadgem2-es_rcp8p5_1950-2099_USA_mon.nc4'.format(dat_dir))

gs      = [3, 4, 5, 6, 7, 8]

da   = da.sel(time=da.time.dt.month.isin(gs)).resample(time='1A').mean()

da.to_netcdf('{0}/tas_bced_1960_1999_hadgem2-es_historical_1950-2099_USA_ann.nc4'.format(dat_dir))


