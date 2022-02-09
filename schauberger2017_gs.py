#!/bin/env python
# from h5py import File as hfile
import pandas as pd
import xarray as xr
import numpy as np
from growingSeasonVar import *


def growingSeasonBins(crop, loc=None, t_shift=0):
    if loc != None:
        top, bottom, left, right = loc
    else:
        top, bottom, left, right = [0, 360, 0, 720]
    dx = bottom - top
    dy = right - left

    prefix = 'schauberger2017'
    netcdf = crop[:3]

    if (crop == 'maize') or (crop == 'soy'):
        # Remember 0 index, so Jan 1st is really 0th day
        plant = np.ma.zeros((dx, dy)) + 59  # March 1st
        har = np.ma.zeros((dx, dy)) + 243  # August 31st
    elif (crop == 'wheat'):
        plant = np.ma.zeros((dx, dy)) + 287  # October 15th
        har = np.ma.zeros((dx, dy)) + 561  # July 15th of next year

    pl = (plant).astype(int)
    pl = np.repeat(pl[np.newaxis, :, :], 730, axis=0)
    ha = (har).astype(int)
    ha = np.repeat(ha[np.newaxis, :, :], 730, axis=0)

    pl = np.where(ha < 365, pl + 365, pl)
    ha = np.where(ha < 365, ha + 365, ha)

    dayids = np.reshape(np.repeat(np.repeat(np.arange(0, 730, 1), dx, axis=np.newaxis), dy, axis=np.newaxis),
                        (730, dx, dy))

    f1 = '/project2/geos39650/ag_data/agmerra/'

    bins = np.arange(0, 43, 3)
    temp_out = np.ma.zeros((30, dx, dy, (bins.size - 1)))
    startyear = 1981  # this is index 0 in the following loop
    leapyear = 1984  # this is the first leap year, following leaps are +4
    leap = 0  # add a day for each leap year cumulatively through the time loop
    for i in range(30):
        print(startyear + i)
        if (i + 1) % 4 == 0:
            leap += 1
            print('leap year')

        tasmini = hfile(f1 + 'tasmin' + '_agmerra_1980-2010.nc4', 'r')['tasmin'][
                  ((i) * 365 + leap):((i) * 365 + 730 + leap), top:bottom, left:right]
        tasmini = np.ma.masked_where(dayids < pl, tasmini)
        tasmini = np.ma.masked_where(dayids > ha, tasmini)
        tasmini = np.ma.masked_where(tasmini == 1e20, tasmini)
        tasmini = np.ma.harden_mask(tasmini)

        tasmaxi = hfile(f1 + 'tasmax' + '_agmerra_1980-2010.nc4', 'r')['tasmax'][
                  ((i) * 365 + leap):((i) * 365 + 730 + leap), top:bottom, left:right]
        tasmaxi = np.ma.masked_where(dayids < pl, tasmaxi)
        tasmaxi = np.ma.masked_where(dayids > ha, tasmaxi)
        tasmaxi = np.ma.masked_where(tasmaxi == 1e20, tasmaxi)
        tasmaxi = np.ma.harden_mask(tasmaxi)

        tasmaxi = np.reshape(np.repeat(tasmaxi, 24, axis=np.newaxis), (730, dx, dy, 24))
        tasmini = np.reshape(np.repeat(tasmini, 24, axis=np.newaxis), (730, dx, dy, 24))
        hrs = np.reshape(np.tile(np.arange(24), (730, dx, dy)), (730, dx, dy, 24))

        cos_hrs = np.cos(hrs * np.pi / 12)
        hrs = 0
        amplitude = (-1) * (tasmaxi - tasmini) / 2
        offset = (tasmaxi + tasmini) / 2
        tasmaxi, tasmini = 0, 0
        t_hrs = (amplitude * cos_hrs) + offset
        amplitude, cos_hrs, offset = 0, 0, 0
        t_hrs = t_hrs - 273.15 + t_shift
        ## BIN CALCULATION (1D)
        t_hrs = np.clip(t_hrs, a_min=t_hrs.min(), a_max=42)

        temp = np.apply_along_axis(lambda a: np.histogram(a, bins=bins)[0], 3, t_hrs)

        t_hrs = 0
        temp_out[i, :, :] = np.sum(temp, axis=0) / 24
        temp = 0
    temp_out[temp_out.mask] = 1e20
    growdir = '/project2/geos39650/ag_data/growseasons/'
    np.save(growdir + '{0}_{1}_tbins_{2}.npy'.format(prefix, netcdf, t_shift), temp_out.data)
    print('binned temperature | saved')
    return temp_out


def growingSeasonXR(region):
    """
    Input arrays should have the correct cropping already, so no need to trim in this function
    :param region: name of region
    :return:
    """

    prefix = 'shau'
    crop = 'maize'

    if crop == 'winter_wheat':
        netcdf = 'wwh'
    elif crop == 'spring_wheat':
        netcdf = 'swh'
    else:
        netcdf = crop[:3]

    if historical:
        hist_ext = 'historical_1950-2004'
        startyear   = 1951 # b/c the pl/har files start with 1951
        endyear     = 2004
    else:
        hist_ext = 'rcp8p5_2005-2099'
        startyear   = 2005
        endyear     = 2099

    # The directory in which the climate data files are located
    f1 = '/project2/geos39650/ag_data/climate_projections/lpjml/'

    tasmin = xr.open_dataarray(f1 + 'tasmin_bced_1960_1999_hadgem2-es_historical_1950-2004_USA.nc4')
    tasmax = xr.open_dataarray(f1 + 'tasmax_bced_1960_1999_hadgem2-es_historical_1950-2004_USA.nc4')
    dx = tasmax.lat.size
    dy = tasmax.lon.size
    startyear = pd.to_datetime(tasmax.time.values[0]).year
    endyear   = pd.to_datetime(tasmax.time.values[-1]).year

    bins = np.arange(0, 43, 3)
    temp_out = np.ma.zeros((endyear - startyear + 1, dx, dy, (bins.size - 1)))
    for i in range(endyear - startyear + 1):
        year = startyear + i
        print(year)

        tasmini = tasmin.sel(time=slice(str(year)+plant_date, str(year)+har_date)).values
        tasmini[np.isnan(tasmini)] = np.NaN
        tasmini = np.ma.harden_mask(tasmini)

        tasmaxi = tasmax.sel(time=slice(str(year)+plant_date, str(year)+har_date)).values
        tasmaxi[np.isnan(tasmaxi)] = np.NaN
        tasmaxi = np.ma.harden_mask(tasmaxi)

        assert tasmaxi.shape == tasmini.shape, "SHAPES DON'T MATCH!"
        dt = tasmaxi.shape[0]

        # Broadcast arrays in 24-hrs and get cosine functions of hourly temperatures
        tasmaxi = np.reshape(np.repeat(tasmaxi, 24, axis=np.newaxis), (dt, dx, dy, 24))
        tasmini = np.reshape(np.repeat(tasmini, 24, axis=np.newaxis), (dt, dx, dy, 24))
        hrs = np.reshape(np.tile(np.arange(24), (dt, dx, dy)), (dt, dx, dy, 24))
        cos_hrs = np.cos(hrs * np.pi / 12)
        hrs = 0

        amplitude = (-1) * (tasmaxi - tasmini) / 2
        offset = (tasmaxi + tasmini) / 2
        tasmaxi, tasmini = 0, 0
        t_hrs = (amplitude * cos_hrs) + offset
        amplitude, cos_hrs, offset = 0, 0, 0
        t_hrs = t_hrs - 273.15
        ## BIN CALCULATION (1D)
        t_hrs = np.clip(t_hrs, a_min=t_hrs.min(), a_max=42)
        temp = np.apply_along_axis(lambda a: np.histogram(a, bins=bins)[0], 3, t_hrs)


        t_hrs = 0
        temp_out[i, :, :] = np.sum(temp, axis=0) / 24
        temp = 0
    temp_out[temp_out.mask] = 1e20
    growdir = '/project2/geos39650/ag_data/growseasons/'
    np.save(growdir + '{0}_{1}_tbins_{2}-{3}_{4}.npy'.format(prefix, netcdf, startyear, endyear, region), temp_out.data)
    print('binned temperature | saved')
    return temp_out


def trueGS_bins(crop, model, t_shift, rainfed=True):
    """
    :param temp_var: 'tbins'
    :param crop:
    :param model:
    :param t_shift:
    :param rainfed:
    :return:
    """
    temp_var = 'tbins'
    netcdf = crop[:3]

    base_path = '/project2/ggcmi/AgMIP.output/{0}/phase2/{1}/A0'.format(model, crop)

    temp_var_path = '/project2/geos39650/ag_data/growseasons/tbins_USA_{1}.nc'.format(temp_var, t_shift)

    temp        = xr.open_dataarray(temp_var_path)
    temp_lat    = temp.lat
    temp_lon    = temp.lon

    temp_out    = xr.DataArray(dims=['time', 'lat', 'lon', 'tbin'],
                               coords={'time':pd.date_range('1981', '2011', freq='1A'),
                                       'lat':temp_lat,
                                       'lon':temp_lon,
                                       'tbin':np.arange(0,14)})

    w = 0 if rainfed else 'inf'
    N = 200 if model != 'CARAIB' else 'NA'
    plant_file  = "{0}/plant-day/{1}_agmerra_fullharm_plant-day_mai_global_annual_1980_2010_C360_T{2}_W{3}_N{4}_A0.nc4".format(base_path, model.lower(), t_shift, w, N)
    plant       = xr.open_dataarray(plant_file, decode_times=False).sel(lat=temp_lat, lon=temp_lon)
    plant       = np.ma.masked_array(plant)
    plant       = np.ma.masked_where(np.isnan(plant), plant)
    plant       = np.ma.harden_mask(plant).astype(int)

    har_file    = "{0}/maty-day/{1}_agmerra_fullharm_maty-day_mai_global_annual_1980_2010_C360_T{2}_W{3}_N{4}_A0.nc4".format(base_path, model.lower(), t_shift, w, N)
    har         = xr.open_dataarray(har_file, decode_times=False).sel(lat=temp_lat, lon=temp_lon)
    har         = np.ma.masked_array(har)
    har         = np.ma.masked_where(np.isnan(har), har)
    har         = np.ma.harden_mask(har).astype(int)

    startyear = 1980
    for i in range(0, 30):
        curryear = startyear + i
        print(curryear, '-', curryear + 1)

        pl = plant[i:i+2].copy()
        pl[1] = pl[1] + 365
        ha = har[i:i+2].copy()
        ha = pl + ha

        temp_i      = temp.sel(time=slice(str(curryear), str(curryear + 1)))
        lenyear1    = temp_i.sel(time=str(curryear)).time.size
        lenyear2    = temp_i.time.size
        dx          = temp_i.lat.size
        dy          = temp_i.lon.size
        temp_i      = temp_i.values

        dayids = np.reshape(
            np.repeat(
                np.repeat(
                    np.repeat(
                        np.arange(1, lenyear2+1, 1), 14, axis=np.newaxis),
                    dx, axis=np.newaxis),
                dy, axis=np.newaxis),
            (lenyear2, dx, dy, 14))

        # Cut harvest dates before the year started
        pl = np.ma.masked_where(ha < lenyear1, pl)
        ha = np.ma.masked_where(ha < lenyear1, ha)
        # Cut harvest dates after the year ended
        pl = np.ma.masked_where(ha > lenyear2, pl)
        ha = np.ma.masked_where(ha > lenyear2, ha)

        # Check if you've got any double-dates
        plsum = np.sum(pl / pl, axis=0)
        plwhere = np.where(plsum > 1)
        hasum = np.sum(ha / ha, axis=0)
        hawhere = np.where(hasum > 1)
        if not ((plsum[plwhere].size == 0) and (hasum[hawhere].size == 0)):
            print('Size', plwhere[0].size)
            for j in range(plwhere[0].size):
                assert plwhere[0][j] == hawhere[0][j], "Double-dates don't match for PL and HA"
                assert plwhere[1][j] == hawhere[1][j], "Double-dates don't match for PL and HA"
                lat_i = plwhere[0][j]
                lon_i = plwhere[1][j]
                pl.mask[0, lat_i, lon_i] = True
                ha.mask[0, lat_i, lon_i] = True

        # Collapse the 2 1-years into a single 2-year 2D array
        pl = np.sum(pl, axis=0)
        ha = np.sum(ha, axis=0)

        # Repeat out to fit daily size
        pl = np.repeat(np.repeat(pl[np.newaxis, :, :], lenyear2, axis=0)[:, :, :, np.newaxis], 14, axis=3)
        ha = np.repeat(np.repeat(ha[np.newaxis, :, :], lenyear2, axis=0)[:, :, :, np.newaxis], 14, axis=3)

        # Assert that same grid cells are masked in PL and HA arrays
        assert np.array_equal(ha.mask, pl.mask), "Masks of dates don't match for PL and HA"

        temp_i = np.ma.masked_where(dayids < pl, temp_i)
        temp_i = np.ma.masked_where(dayids > ha, temp_i)
        temp_i = np.ma.masked_where(temp_i >= 1e20, temp_i)
        temp_i = np.ma.harden_mask(temp_i)

        temp_out.loc[dict(time=str(curryear + 1))] = np.sum(temp_i, axis=0)

    save_dir = '/project2/geos39650/ag_data/true_gs/'
    temp_out.to_netcdf('{0}{1}_{2}_TrueGS_{3}_T{4}_W{5}.nc'.format(save_dir, model.lower(), netcdf, temp_var, t_shift, w))
    return True


def RCPgs():
    crop = 'maize'
    netcdf = crop[:3]

    # # get GDD and HDD
    # growingSeasonXR('USA')


    # get pr values
    var = 'pr'
    aggfxn = lambda x: np.ma.sum(x, axis=0)
    save_name = '{0}_{1}_pr_{2}-{3}_{4}'.format('shau', netcdf, 2031, 2099, 'USA')
    growingSeasonVarXR(var, aggfxn, '-03-01', '-08-31', save_name)



def fixedGSlpjmlRCP(temp_var, historical):
    """
    :param temp_var: can be TBINS or PR
    :param crop:
    :param model:
    :param t_shift:
    :return:
    """
    crop = 'maize'
    netcdf = crop[:3]
    model = 'LPJmL'
    plant_date = '-03-01'
    har_date = '-08-31'

    if historical:
        hist_ext = 'historical_1950-2004'
        startyear   = 1951 # b/c the pl/har files start with 1951
        endyear     = 2004
    else:
        hist_ext = 'rcp8p5_2005-2099'
        startyear   = 2005
        endyear     = 2099

    if temp_var == 'pr':
        temp_var_path = '/project2/geos39650/ag_data/climate_projections/lpjml/pr_bced_1960_1999_hadgem2-es_{0}_USA.nc4'.format(hist_ext)
        temp = xr.open_dataarray(temp_var_path)
    elif temp_var == 'tbins':
        temp_var_path = '/project2/geos39650/ag_data/growseasons/rcp_{0}_{1}_USA.nc4'.format(temp_var, hist_ext)
        temp = xr.open_dataarray(temp_var_path)

    if historical:
        temp = temp.sel(time=slice('1951', '2004'))

    temp_lat     = temp.lat
    temp_lon     = temp.lon

    if temp_var == 'pr':
        temp_out = xr.DataArray(dims=['time', 'lat', 'lon'],
                                coords={'time': pd.date_range(str(startyear + 1), str(endyear + 1), freq='1A'),
                                        'lat': temp_lat,
                                        'lon': temp_lon})
    elif temp_var == 'tbins':
        temp_out = xr.DataArray(dims=['time', 'lat', 'lon', 'tbin'],
                                coords={'time': pd.date_range(str(startyear + 1), str(endyear + 1), freq='1A'),
                                        'lat': temp_lat,
                                        'lon': temp_lon,
                                        'tbin': np.arange(14)})

    dy = endyear - startyear
    for i in range(1, dy + 1):
        year = startyear + i
        print(year)

        temp_i = temp.sel(time=slice(str(year) + plant_date, str(year) + har_date)).values
        temp_i = np.ma.masked_where(np.isnan(temp_i), temp_i)
        temp_out.loc[dict(time=str(year))] = np.ma.sum(temp_i, axis=0)

    save_dir = '/project2/geos39650/ag_data/growseasons/'
    temp_out.to_netcdf('{0}{1}_{2}_RCP_shau_{3}_{4}_fixed.nc'.format(save_dir, model.lower(), netcdf, hist_ext, temp_var))
    return True



if __name__ == '__main__':
    one = fixedGSlpjmlRCP('tbins', historical=True)
    two = fixedGSlpjmlRCP('tbins', historical=False)
    three = fixedGSlpjmlRCP('pr', historical=True)
    four = fixedGSlpjmlRCP('pr', historical=False)

