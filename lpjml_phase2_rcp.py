import xarray as xr
import numpy as np
import pandas as pd

##########################
###### FUNCTIONS #########
##########################

def compileOutput(output_var, coords=[(48.75, 36.25), (-103.75, -80.25)]):
    output_dir = '/project2/ggcmi/AgMIP.output/LPJmL/phase2/fasttrack/HadGEM2-ES/rcp8p5_v2_201120/rcp8p5_v2_201120/mai/c360_n200'
    fname = '{0}/lpjml_hadgem2-es_fullharm_{1}_mai_global_annual_1951_2099_noirr2.nc4'.format(output_dir,
                                                                                              output_var)
    if coords != None:
        lat_N, lat_S = coords[0]
        lon_W, lon_E = coords[1]
        da = xr.open_dataarray(fname, decode_times=False).sel(lat=slice(lat_N, lat_S), lon=slice(lon_W, lon_E))
    else:
        da = xr.open_dataarray(fname, decode_times=False)
    da['time'] = pd.date_range('1951', '2100', freq='1A')
    da = da.sel(time=slice('1971', '2099'))
    savedir = '/project2/geos39650/ag_data/ggcmi_phase2'
    da.to_netcdf('{0}/lpjml_hadgem2-es_fullharm_c360_n200_{1}_mai_global_annual_1971_2099_noirr2.nc'.format(savedir, output_var))
    return True

def loadDailyData(clim_var, coords=[(48.75, 36.25), (-103.75, -80.25)]):
    if (clim_var == 'tasmin') or (clim_var == 'tasmax') or (clim_var == 'tas'):
        clim_ext = clim_var + '_v2'
    else:
        clim_ext = clim_var + '_v3'
    clim_dir = '/project2/ggcmi/ISIMIP_fasttrack/HadGEM2-ES'
    years = ['1971-1980',
             '1981-1990',
             '1991-2000',
             '2001-2004',
             '2005-2010',
             '2011-2020',
             '2021-2030',
             '2031-2040',
             '2041-2050',
             '2051-2060',
             '2061-2070',
             '2071-2080',
             '2081-2090',
             '2091-2099']
    scen = ['historical',
            'historical',
            'historical',
            'historical',
            'rcp8p5',
            'rcp8p5',
            'rcp8p5',
            'rcp8p5',
            'rcp8p5',
            'rcp8p5',
            'rcp8p5',
            'rcp8p5',
            'rcp8p5',
            'rcp8p5']
    scen = dict(zip(years, scen))

    file_names = ['{0}/{1}/{2}/{3}_bced_1960_1999_hadgem2-es_{1}_{4}.nc4'.format(clim_dir,
                                                                                 scen[y],
                                                                                 clim_ext,
                                                                                 clim_var,
                                                                                 y) for y in years]
    if coords != None:
        lat_N, lat_S = coords[0]
        lon_W, lon_E = coords[1]
        da = xr.concat(
            [xr.open_dataarray(fname).sel(lat=slice(lat_N, lat_S), lon=slice(lon_W, lon_E)) for fname in file_names],
            dim='time')
    else:
        da = xr.concat([xr.open_dataarray(fname) for fname in file_names], dim='time')
    return da


def loadPR():
    pr = loadDailyData('pr')
    save_dir = '/project2/geos39650/ag_data/climate_projections/lpjml'
    pr.to_netcdf('{0}/lpjml_rcp8p5_{1}_MDW_daily_1971_2099.nc'.format(save_dir, 'pr'))
    return True


def loadGDDandHDD():
    tasmin = loadDailyData('tasmin')
    tasmax = loadDailyData('tasmax')

    gdd = xr.DataArray(np.zeros(tasmin.shape),
                       coords={'time': tasmin.time, 'lat': tasmin.lat, 'lon': tasmin.lon},
                       dims=['time', 'lat', 'lon'])
    hdd = xr.DataArray(np.zeros(tasmin.shape),
                       coords={'time': tasmin.time, 'lat': tasmin.lat, 'lon': tasmin.lon},
                       dims=['time', 'lat', 'lon'])

    for y in np.sort(np.unique(tasmin['time.year'].values)):
        print(y)
        tasmini = tasmin.sel(time=str(y))
        tasmaxi = tasmax.sel(time=str(y))

        dt, dx, dy = tasmini.shape

        tasmini_vals = np.ma.masked_array(np.reshape(np.repeat(tasmini.values, 24, axis=np.newaxis), (dt, dx, dy, 24)))
        tasmini_vals = np.ma.masked_where(np.isnan(tasmini_vals), tasmini_vals)
        tasmini_vals = np.ma.harden_mask(tasmini_vals)

        tasmaxi_vals = np.ma.masked_array(np.reshape(np.repeat(tasmaxi.values, 24, axis=np.newaxis), (dt, dx, dy, 24)))
        tasmaxi_vals = np.ma.masked_where(np.isnan(tasmaxi_vals), tasmaxi_vals)
        tasmaxi_vals = np.ma.harden_mask(tasmaxi_vals)

        hrs = np.reshape(np.tile(np.arange(24), (dt, dx, dy)), (dt, dx, dy, 24))

        cos_hrs = np.cos(hrs * np.pi / 12)
        amplitude = (-1) * (tasmaxi_vals - tasmini_vals) / 2
        offset = (tasmaxi_vals + tasmini_vals) / 2
        t_hrs = (amplitude * cos_hrs) + offset
        t_hrs = t_hrs - 273.15
        gdd_hrs = t_hrs.copy()
        hdd_hrs = t_hrs.copy()

        t_high = 29
        t_low = 10
        ## GDD CALCULATION
        gdd_hrs[gdd_hrs > t_high] = t_high
        gdd_hrs = gdd_hrs - t_low
        gdd_hrs[gdd_hrs < 0] = 0
        gdd.loc[dict(time=tasmini.time)] = np.sum(gdd_hrs * (1 / 24), axis=3)

        ## HDD CALCULATION
        hdd_hrs = hdd_hrs - t_high
        hdd_hrs[hdd_hrs < 0] = 0
        hdd.loc[dict(time=tasmini.time)] = np.sum(hdd_hrs * (1 / 24), axis=3)

    save_dir = '/project2/geos39650/ag_data/climate_projections/lpjml'
    gdd.to_netcdf('{0}/lpjml_rcp8p5_{1}_MDW_daily_1971_2099.nc'.format(save_dir, 'gdd'))
    hdd.to_netcdf('{0}/lpjml_rcp8p5_{1}_MDW_daily_1971_2099.nc'.format(save_dir, 'hdd'))
    return True


def fixedGS_rcp(clim_var):
    """
    A function to compile the fixed-season climate variables
    :param clim_var: name of wanted variable
    :return: True after making fixed-season file
    """
    plant_date = '-03-01'
    har_date = '-08-28'

    clim_dir = '/project2/geos39650/ag_data/climate_projections/lpjml'
    clim_var_path = '{0}/lpjml_rcp8p5_{1}_MDW_daily_1971_2099.nc'.format(clim_dir, clim_var)
    clim_arr = xr.open_dataarray(clim_var_path)

    years = clim_arr['time.year'].resample(time='1A').max().values

    out  = xr.DataArray(dims=['time', 'lat', 'lon'],
                        coords={'time':pd.date_range(str(years[0]), str(years[-1]+1), freq='1A'),
                                'lat':clim_arr.lat,
                                'lon':clim_arr.lon})

    for year in years:
        print(year)

        clim_year = clim_arr.sel(time=slice(str(year) + plant_date, str(year) + har_date)).values
        clim_year = np.ma.masked_where(np.isnan(clim_year), clim_year)
        out.loc[dict(time=str(year))] = np.ma.sum(clim_year, axis=0)

    save_dir = '/project2/geos39650/ag_data/growseasons/lpjml'
    out.to_netcdf('{0}/lpjml_rcp8p5_{1}_MDW_fixedGS_1971_2099.nc'.format(save_dir, clim_var))
    return True

def trueGS_rcp(clim_var):
    """
    Get the true-season values for GDD, HDD, Pr
    :param clim_var: variable of climate to calculate
    :return:
    """
    startyear = 1971
    endyear   = 2099

    clim_dir = '/project2/geos39650/ag_data/climate_projections/lpjml'
    clim_var_path = '{0}/lpjml_rcp8p5_{1}_MDW_daily_1971_2099.nc'.format(clim_dir, clim_var)
    clim_arr = xr.open_dataarray(clim_var_path)
    extra_year = xr.DataArray(dims=['time', 'lat', 'lon'],
                            coords={'time': pd.date_range(str(endyear+1)+'-01-01', str(endyear+1)+'-12-31', freq='1D'),
                                    'lat': clim_arr.lat,
                                    'lon': clim_arr.lon})
    clim_arr = xr.concat([clim_arr, extra_year], dim='time')
    out = xr.DataArray(dims=['time', 'lat', 'lon'],
                       coords={'time':pd.date_range(str(startyear), str(endyear+1), freq='1A'),
                               'lat':clim_arr.lat,
                               'lon':clim_arr.lon})

    PlHa_dir = '/project2/geos39650/ag_data/ggcmi_phase2'
    plant_file  = '{0}/lpjml_hadgem2-es_fullharm_c360_n200_plantday_mai_global_annual_1971_2099_noirr2.nc4'.format(PlHa_dir)
    plant       = xr.open_dataarray(plant_file).values
    plant       = np.append(plant, np.empty(plant[0].shape)[np.newaxis,:, :]*np.NaN, axis=0)
    plant       = np.ma.masked_array(plant)
    plant       = np.ma.masked_where(np.isnan(plant), plant)
    plant       = np.ma.harden_mask(plant).astype(int)

    har_file    = '{0}/lpjml_hadgem2-es_fullharm_c360_n200_matyday_mai_global_annual_1971_2099_noirr2.nc4'.format(PlHa_dir)
    har         = xr.open_dataarray(har_file).values
    har         = np.append(har, np.empty(har[0].shape)[np.newaxis, :, :]*np.NaN, axis=0)
    har         = np.ma.masked_array(har)
    har         = np.ma.masked_where(np.isnan(har), har)
    har         = np.ma.harden_mask(har).astype(int)

    for i in range(0, endyear - startyear + 1):
        curryear = startyear + i
        print(curryear, '-', curryear + 1)
        clim_year       = clim_arr.sel(time=slice(str(curryear), str(curryear + 1)))
        lenyear1        = clim_year.sel(time=str(curryear)).time.size
        lenyear2        = clim_year.time.size
        dx              = clim_year.lat.size
        dy              = clim_year.lon.size
        clim_year       = clim_year.values

        pl = plant[i:i + 2].copy()
        pl = np.ma.masked_where(np.isnan(pl), pl)
        # Mask PLANTING dates WHERE ZERO (never planted)
        pl = np.ma.masked_where(pl == 0, pl)
        pl[1] = pl[1] + lenyear1

        ha = har[i:i + 2].copy()
        ha = np.ma.masked_where(np.isnan(ha), ha)
        # Mask HARVEST dates WHERE ZERO (failure)
        ha = np.ma.masked_where(ha == 0, ha)
        ha = pl + ha

        dayids = np.reshape(np.repeat(np.repeat(np.arange(1, lenyear2+1, 1), dx, axis=np.newaxis), dy, axis=np.newaxis),
                            (lenyear2, dx, dy))

        # Mask PLANTING dates AFTER Year-1 ended
        pl = np.ma.masked_where(pl > lenyear1, pl)
        ha = np.ma.masked_where(pl > lenyear1, ha)

        # Mask HARVEST dates AFTER Year-2 ended
        pl = np.ma.masked_where(ha > lenyear2, pl)
        ha = np.ma.masked_where(ha > lenyear2, ha)

        if i == dy-1: # If it's the last year of calculation
            # Cut HARVEST dates AFTER the last REAL year (Year-1)
            pl = np.ma.masked_where(ha > lenyear1, pl)
            ha = np.ma.masked_where(ha > lenyear1, ha)

        # Check if you've got any double-dates
        plsum = np.ma.sum(pl / pl, axis=0)
        plwhere = np.where(plsum > 1)
        hasum = np.ma.sum(ha / ha, axis=0)
        hawhere = np.where(hasum > 1)
        if not ((plsum[plwhere].size == 0) and (hasum[hawhere].size == 0)):
            print('Size', plwhere[0].size)
            for j in range(plwhere[0].size):
                assert plwhere[0][j] == hawhere[0][j], "Double-dates don't match for PL and HA"
                assert plwhere[1][j] == hawhere[1][j], "Double-dates don't match for PL and HA"
                lat_i = plwhere[0][j]
                lon_i = plwhere[1][j]
                print("[{0}, {1}]".format(clim_arr.lat.values[lat_i], clim_arr.lon.values[lon_i]))
                pl.mask[1, lat_i, lon_i] = True     # Mask the earlier harvest, only show the latter
                ha.mask[1, lat_i, lon_i] = True     # Mask the earlier harvest, only show the latter


        # Collapse the 2 1-years into a single 2-year 2D array
        pl = np.sum(pl, axis=0)
        ha = np.sum(ha, axis=0)

        # Repeat out to fit daily size
        pl = np.repeat(pl[np.newaxis, :, :], lenyear2, axis=0)
        ha = np.repeat(ha[np.newaxis, :, :], lenyear2, axis=0)

        # Assert that same grid cells are masked in PL and HA arrays
        assert np.array_equal(ha.mask, pl.mask), "Masks of dates don't match for PL and HA"

        clim_year = np.ma.masked_where(dayids < pl, clim_year)
        clim_year = np.ma.masked_where(dayids > ha, clim_year)
        clim_year = np.ma.masked_where(clim_year >= 1e20, clim_year)
        clim_year = np.ma.harden_mask(clim_year)

        out.loc[dict(time=str(curryear))] = np.ma.sum(clim_year, axis=0)

    save_dir = '/project2/geos39650/ag_data/growseasons/lpjml'
    out.to_netcdf('{0}/lpjml_rcp8p5_{1}_MDW_trueGS_1971_2099.nc'.format(save_dir, clim_var))
    return True


##########################
###### EXECUTION #########
##########################

def getDailyValues(): # Ran on Nov. 24
    print('loading PR')
    loadPR()
    print('LOADED PR')
    print('loading TEMP')
    loadGDDandHDD()
    print('LOADED TEMP')
    return True

def compile_GGCMI_files(): # Ran on Nov. 24
    output_vars = ['matyday', 'plantday', 'yield']
    for output_var in output_vars:
        print('VARIABLE: {0}'.format(output_var))
        print('\n')
        output_loaded = compileOutput(output_var)
        print('COMPLETE: {0}'.format(output_loaded))

def get_GSmean_tas(): # Ran on Nov. 24
    print('loading TAS')
    tas = loadDailyData('tas')
    tas = tas.sel(time=tas.time.dt.month.isin([3,4,5,6,7,8])).resample(time="1A").mean()
    save_dir = '/project2/geos39650/ag_data/climate_projections/lpjml'
    tas.to_netcdf('{0}/lpjml_rcp8p5_{1}_MDW_GSmean_1971_2099.nc'.format(save_dir, 'tas'))
    print('LOADED TAS')
    return True

def getSeasonSums(): # Ran on Nov. 24
    clim_vars = ['gdd', 'hdd', 'pr']
    for clim_var in clim_vars:
        print('VARIABLE: {0}'.format(clim_var))
        clim_loaded = fixedGS_rcp(clim_var)
        print('COMPLETE: {0}'.format(clim_loaded))
        print('\n')

def get_trueGS_Sums(): # Ran on Nov. 24
    clim_vars = ['gdd', 'hdd', 'pr']
    for clim_var in clim_vars:
        print('VARIABLE: {0}'.format(clim_var))
        clim_loaded = trueGS_rcp(clim_var)
        print('COMPLETE: {0}'.format(clim_loaded))
        print('\n')


#####################
###### MAIN #########
#####################
if __name__ == '__main__':
    get_trueGS_Sums()

