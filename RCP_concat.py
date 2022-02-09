import xarray as xr

def concatFiles(var, file_path, file_list, outname):
    top = 49.3457868  # north lat
    left = -124.7844079  # west long
    right = -66.9513812  # east long
    bottom = 24.7433195  # south lat

    file_names = [file_path + var + name for name in file_list]
    print('NAMES DONE.')

    new_ds = xr.concat([xr.open_dataset(name)[var].sel(lat=slice(top,bottom), lon=slice(left,right)) for name in file_names], dim='time').sortby('time')
    new_ds = xr.concat([xr.open_dataset(name)[var].sel(lat=42.25, lon=-93.75, method='nearest') for name in file_names], dim='time').sortby('time')

    print('CONCAT DONE.')

    new_ds.to_netcdf(outname)
    print(outname)
    print('COMPLETE.')

vars = ['pr', 'tasmax', 'tasmin']
exts = ['v3', 'v2', 'v2']
for var, ext in zip(vars, exts):
    file_path = '/project2/ggcmi/ISIMIP_fasttrack/HadGEM2-ES/rcp8p5/{0}_{1}/'.format(var, ext)
    file_list = ['_bced_1960_1999_hadgem2-es_rcp8p5_2005-2010.nc4',
                 '_bced_1960_1999_hadgem2-es_rcp8p5_2011-2020.nc4',
                 '_bced_1960_1999_hadgem2-es_rcp8p5_2021-2030.nc4',
                 '_bced_1960_1999_hadgem2-es_rcp8p5_2031-2040.nc4',
                 '_bced_1960_1999_hadgem2-es_rcp8p5_2041-2050.nc4',
                 '_bced_1960_1999_hadgem2-es_rcp8p5_2051-2060.nc4',
                 '_bced_1960_1999_hadgem2-es_rcp8p5_2061-2070.nc4',
                 '_bced_1960_1999_hadgem2-es_rcp8p5_2071-2080.nc4',
                 '_bced_1960_1999_hadgem2-es_rcp8p5_2081-2090.nc4',
                 '_bced_1960_1999_hadgem2-es_rcp8p5_2091-2099.nc4']


    outname   = '/project2/geos39650/ag_data/climate_projections/lpjml/{0}_bced_1960_1999_hadgem2-es_rcp8p5_2005-2099_USA.nc4'.format(var)
    concatFiles(var, file_path, file_list, outname)
