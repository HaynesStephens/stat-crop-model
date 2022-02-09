import xarray as xr
import numpy as np

def getHistTemps():
    top = 49.3457868  # north lat
    left = -124.7844079  # west long
    right = -66.9513812  # east long
    bottom = 24.7433195  # south lat
    var = 'tas'

    histdir = '/project2/ggcmi/ISIMIP_fasttrack/HadGEM2-ES/historical/tas_v2/'
    ncfiles = ['tas_bced_1960_1999_hadgem2-es_historical_1981-1990.nc4',
               'tas_bced_1960_1999_hadgem2-es_historical_1991-2000.nc4',
               'tas_bced_1960_1999_hadgem2-es_historical_2001-2004.nc4']
    file_names = [histdir + file_i for file_i in ncfiles]
    ds       = xr.concat([xr.open_dataset(name)[var].sel(lat=slice(top,bottom), lon=slice(left,right), time=slice('1981','2004')) for name in file_names], dim='time').sortby('time')
    ds       = ds.sel(time=ds.time.dt.month.isin([3,4,5,6,7,8]))
    ds_r     = ds.resample(time="1A").mean()

    gs_T = np.array([288.78064, 289.10275, 287.9848 , 288.1931 , 288.31223, 289.5408 ,
                     289.15308, 289.62558, 289.59143, 289.69443, 289.41837, 288.66492,
                     288.95792, 289.6768 , 289.69464, 289.852  , 290.7099 , 290.317  ,
                     289.7323 , 289.72528, 289.31158, 290.3329 , 290.50586, 290.26602])

    gs_meanT = 289.46436
    return gs_T, gs_meanT


top = 49.3457868  # north lat
left = -124.7844079  # west long
right = -66.9513812  # east long
bottom = 24.7433195  # south lat
var = 'tas'

histdir = '/project2/ggcmi/ISIMIP_fasttrack/HadGEM2-ES/rcp8p5/tas_v2/'
ncfiles = ['tas_bced_1960_1999_hadgem2-es_rcp8p5_2031-2040.nc4',
           'tas_bced_1960_1999_hadgem2-es_rcp8p5_2041-2050.nc4',
           'tas_bced_1960_1999_hadgem2-es_rcp8p5_2051-2060.nc4',
           'tas_bced_1960_1999_hadgem2-es_rcp8p5_2061-2070.nc4',
           'tas_bced_1960_1999_hadgem2-es_rcp8p5_2071-2080.nc4',
           'tas_bced_1960_1999_hadgem2-es_rcp8p5_2081-2090.nc4',
           'tas_bced_1960_1999_hadgem2-es_rcp8p5_2091-2099.nc4']

file_names = [histdir + file_i for file_i in ncfiles]
ds       = xr.concat([xr.open_dataset(name)[var].sel(lat=slice(top,bottom), lon=slice(left,right), time=slice('2031','2099')) for name in file_names], dim='time').sortby('time')
ds       = ds.sel(time=ds.time.dt.month.isin([3,4,5,6,7,8]))
ds_r     = ds.resample(time="1A").mean()

gs_t = np.array([291.06763, 291.15308, 292.14502, 292.1896 , 291.7429 ,
                 291.84903, 291.29645, 291.835  , 291.3287 , 291.47433,
                 293.0296 , 291.87888, 292.6714 , 292.30014, 292.77597,
                 292.65546, 292.65442, 292.2926 , 292.94824, 292.02966,
                 293.35138, 292.3011 , 293.32336, 293.01184, 293.32288,
                 293.58313, 293.1122 , 293.7403 , 292.8281 , 294.38318,
                 294.1053 , 294.54614, 293.84238, 294.59354, 294.05063,
                 293.8092 , 293.8909 , 294.36002, 293.95   , 294.6826 ,
                 294.69427, 295.05093, 294.84274, 294.59   , 294.37216,
                 295.4597 , 294.46188, 295.8914 , 294.8636 , 295.3261 ,
                 295.43048, 295.62457, 295.72763, 295.65253, 295.23584,
                 295.22678, 295.3174 , 295.3113 , 295.73615, 296.2149 ,
                 295.49878, 295.99677, 296.54498, 296.15784, 296.65674,
                 296.58118, 295.40414, 296.66678, 295.82883])
