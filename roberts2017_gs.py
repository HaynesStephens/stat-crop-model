import numpy as np
# from h5py import File as hfile
import pandas as pd
import xarray as xr
from growingSeasonVar import *

def growingSeasonGDD_region(loc, region, t_shift=0):
    top, bottom, left, right = loc
    dx = bottom - top
    dy = right - left
    print(dx, dy)

    prefix = 'roberts2017'
    crop = 'maize'          # Roberts only does maize, any other crop would require a new growing season definition

    if crop == 'winter_wheat':
        netcdf = 'wwh'
    elif crop == 'spring_wheat':
        netcdf = 'swh'
    else:
        netcdf = crop[:3]

    plant = np.ma.zeros((dx, dy)) + 59  # March 1st
    har = np.ma.zeros((dx, dy)) + 180   # 180 days growing season
    ha = plant + har                    # add the 180 days to get the harvesting d.o.y.

    pl = (plant).astype(int)
    pl = np.repeat(pl[np.newaxis, :, :], 730, axis=0)
    ha = (ha).astype(int)
    ha = np.repeat(ha[np.newaxis, :, :], 730, axis=0)

    pl = np.where(ha < 365, pl + 365, pl)
    ha = np.where(ha < 365, ha + 365, ha)

    dayids = np.reshape(np.repeat(np.repeat(np.arange(0, 730, 1), dx, axis=np.newaxis), dy, axis=np.newaxis),
                        (730, dx, dy))

    # The directory in which the AgMERRA data files are located
    f1 = '/project2/geos39650/ag_data/agmerra/'

    Gdd = np.ma.zeros((30, dx, dy))
    Hdd = np.ma.zeros((30, dx, dy))
    startyear = 1981  # this is index 0 in the following loop
    leapyear = 1984  # this is the first leap year, following leaps are +4
    leap = 0  # add a day for each leap year cumulatively through the time loop
    for i in range(30):
        print(startyear + i)
        if (i + 1) % 4 == 0:
            leap += 1
            print('leap year')

        tasmini = hfile(f1 + 'tasmin' + '_agmerra_1980-2010.nc4', 'r')['tasmin'][((i) * 365 + leap):((i) * 365 + 730 + leap), top:bottom, left:right]
        tasmini = np.ma.masked_where(dayids < pl, tasmini)
        tasmini = np.ma.masked_where(dayids > ha, tasmini)
        tasmini = np.ma.masked_where(tasmini >= 1e20, tasmini)
        tasmini = np.ma.harden_mask(tasmini)

        tasmaxi = hfile(f1 + 'tasmax' + '_agmerra_1980-2010.nc4', 'r')['tasmax'][((i) * 365 + leap):((i) * 365 + 730 + leap), top:bottom, left:right]
        tasmaxi = np.ma.masked_where(dayids < pl, tasmaxi)
        tasmaxi = np.ma.masked_where(dayids > ha, tasmaxi)
        tasmaxi = np.ma.masked_where(tasmaxi >= 1e20, tasmaxi)
        tasmaxi = np.ma.harden_mask(tasmaxi)

        # Broadcast arrays in 24-hrs and get cosine functions of hourly temperatures
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
        gdd_hrs = t_hrs.copy()
        hdd_hrs = t_hrs.copy()
        t_hrs = 0

        t_high = 29
        t_low = 10
        ## GDD CALCULATION
        gdd_hrs[gdd_hrs > t_high] = t_high
        gdd_hrs = gdd_hrs - t_low
        gdd_hrs[gdd_hrs < 0] = 0
        Gdd[i, :, :] = np.sum(np.sum(gdd_hrs * (1 / 24), axis=3), axis=0)
        gdd_hrs = 0

        ## HDD CALCULATION
        hdd_hrs = hdd_hrs - t_high
        hdd_hrs[hdd_hrs < 0] = 0
        Hdd[i, :, :] = np.sum(np.sum(hdd_hrs * (1 / 24), axis=3), axis=0)
        hdd_hrs = 0

    Gdd[Gdd.mask] = 1e20
    growdir = '/project2/geos39650/ag_data/growseasons/'
    np.save(growdir + '{0}_{1}_{2}_gdd_{3}.npy'.format(prefix, region, netcdf, t_shift), Gdd.data)
    print('gdd | saved')
    Hdd[Hdd.mask] = 1e20
    np.save(growdir + '{0}_{1}_{2}_hdd_{3}.npy'.format(prefix, region, netcdf, t_shift), Hdd.data)
    print('hdd | saved')
    return True


def midwestGS():
    crop = 'maize'
    netcdf = crop[:3]

    # get GDD and HDD
    growingSeasonGDD_region(loc = [82, 108, 152, 200], region='midwest')

    # get pr values
    loc = [82, 108, 152, 200]
    var = 'pr'
    plant = np.ma.zeros((360, 720)) + 59  # March 1st
    har = np.ma.zeros((360, 720)) + 180  # 180 days growing season
    ha = plant + har  # add the 180 days to get the harvesting d.o.y.
    pl = (plant).astype(int)
    ha = (ha).astype(int)
    aggfxn = lambda x: np.ma.sum(x, axis=0)
    growingSeasonVar(var, aggfxn, pl, ha, '{0}_{1}_{2}'.format('roberts2017', netcdf, var), loc=loc)


def usaGS():
    Tshifts = [-1, 1, 2, 3, 4, 6]
    crop = 'maize'
    netcdf = crop[:3]
    loc  = [81, 131, 110, 226]

    # get GDD and HDD
    for t_shift in Tshifts:
        growingSeasonGDD_region(loc = loc, region='USA', t_shift = t_shift)

    # # get pr values
    # var = 'pr'
    # plant = np.ma.zeros((360, 720)) + 59  # March 1st
    # har = np.ma.zeros((360, 720)) + 180  # 180 days growing season
    # ha = plant + har  # add the 180 days to get the harvesting d.o.y.
    # pl = (plant).astype(int)
    # ha = (ha).astype(int)
    # aggfxn = lambda x: np.ma.sum(x, axis=0)
    # growingSeasonVar(var, aggfxn, pl, ha, '{0}_{1}_{2}'.format('roberts2017_USA', netcdf, var), loc=loc)


def growingSeasonGDD_global():
    loc = [0, 360, 0, 720]
    top, bottom, left, right = loc
    dx = bottom - top
    dy = right - left
    print(dx, dy)

    prefix = 'roberts2017'
    crop = 'maize'          # Roberts only does maize, any other crop would require a new growing season definition

    if crop == 'winter_wheat':
        netcdf = 'wwh'
    elif crop == 'spring_wheat':
        netcdf = 'swh'
    else:
        netcdf = crop[:3]

    plant = np.ma.zeros((dx, dy)) + 59  # March 1st
    har = np.ma.zeros((dx, dy)) + 180   # 180 days growing season
    ha = plant + har                    # add the 180 days to get the harvesting d.o.y.

    pl = (plant).astype(int)
    pl = np.repeat(pl[np.newaxis, :, :], 730, axis=0)
    ha = (ha).astype(int)
    ha = np.repeat(ha[np.newaxis, :, :], 730, axis=0)

    pl = np.where(ha < 365, pl + 365, pl)
    ha = np.where(ha < 365, ha + 365, ha)

    # FIXXXXXX
    # n_xcuts = 4
    # x_cuts = np.arange(0, 360, 90)
    # dx = 90
    # n_ycuts = 4
    # y_cuts = np.arange(0, 720, 180)
    # dy = 180

    dayids = np.reshape(np.repeat(np.repeat(np.arange(0, 730, 1), dx, axis=np.newaxis), dy, axis=np.newaxis),
                        (730, dx, dy))

    # The directory in which the AgMERRA data files are located
    f1 = '/project2/geos39650/ag_data/agmerra/'

    Gdd = np.ma.zeros((30, dx, dy))
    Hdd = np.ma.zeros((30, dx, dy))
    startyear = 1981  # this is index 0 in the following loop
    leapyear = 1984  # this is the first leap year, following leaps are +4
    leap = 0  # add a day for each leap year cumulatively through the time loop
    for i in range(30):
        print(startyear + i)
        if (i + 1) % 4 == 0:
            leap += 1
            print('leap year')

        for j in range(n_xcuts):
            xcutj = x_cuts[j]
            xendj = xcutj + dx

            for k in range(n_ycuts):
                ycutk = y_cuts[k]
                yendk = ycutk + dy
                print('X | {0}:{1}'.format(xcutj, xendj))
                print('Y | {0}:{1}'.format(ycutk, yendk))
                print('\n')

                plj = pl[:, xcutj:xendj, ycutk:yendk]
                haj = ha[:, xcutj:xendj, ycutk:yendk]

                tasmini = hfile(f1 + 'tasmin' + '_agmerra_1980-2010.nc4', 'r')['tasmin'][((i) * 365 + leap):((i) * 365 + 730 + leap), xcutj:xendj, ycutk:yendk]
                tasmini = np.ma.masked_where(dayids < plj, tasmini)
                tasmini = np.ma.masked_where(dayids > haj, tasmini)
                tasmini = np.ma.masked_where(tasmini >= 1e20, tasmini)
                tasmini = np.ma.harden_mask(tasmini)

                tasmaxi = hfile(f1 + 'tasmax' + '_agmerra_1980-2010.nc4', 'r')['tasmax'][((i) * 365 + leap):((i) * 365 + 730 + leap), xcutj:xendj, ycutk:yendk]
                tasmaxi = np.ma.masked_where(dayids < plj, tasmaxi)
                tasmaxi = np.ma.masked_where(dayids > haj, tasmaxi)
                tasmaxi = np.ma.masked_where(tasmaxi >= 1e20, tasmaxi)
                tasmaxi = np.ma.harden_mask(tasmaxi)

                # Broadcast arrays in 24-hrs and get cosine functions of hourly temperatures
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

                t_hrs = t_hrs - 273.15
                gdd_hrs = t_hrs.copy()
                hdd_hrs = t_hrs.copy()
                t_hrs = 0

                t_high = 29
                t_low = 10
                ## GDD CALCULATION
                gdd_hrs[gdd_hrs > t_high] = t_high
                gdd_hrs = gdd_hrs - t_low
                gdd_hrs[gdd_hrs < 0] = 0
                Gdd[i, xcutj:xendj, ycutk:yendk] = np.sum(np.sum(gdd_hrs * (1 / 24), axis=3), axis=0)
                gdd_hrs = 0

                ## HDD CALCULATION
                hdd_hrs = hdd_hrs - t_high
                hdd_hrs[hdd_hrs < 0] = 0
                Hdd[i, xcutj:xendj, ycutk:yendk] = np.sum(np.sum(hdd_hrs * (1 / 24), axis=3), axis=0)
                hdd_hrs = 0

    Gdd[Gdd.mask] = 1e20
    growdir = '/project2/geos39650/ag_data/growseasons/'
    np.save(growdir + '{0}_{1}_gdd_global.npy'.format(prefix, netcdf), Gdd.data)
    print('gdd | saved')
    Hdd[Hdd.mask] = 1e20
    np.save(growdir + '{0}_{1}_hdd_global.npy'.format(prefix, netcdf), Hdd.data)
    print('hdd | saved')
    return True


def globalGS():
    crop = 'maize'
    netcdf = crop[:3]

    # get GDD and HDD
    growingSeasonGDD_global()

    # get pr values
    var = 'pr'
    plant = np.ma.zeros((360, 720)) + 59  # March 1st
    har = np.ma.zeros((360, 720)) + 180  # 180 days growing season
    ha = plant + har  # add the 180 days to get the harvesting d.o.y.
    pl = (plant).astype(int)
    ha = (ha).astype(int)
    aggfxn = lambda x: np.ma.sum(x, axis=0)
    growingSeasonVar(var, aggfxn, pl, ha, '{0}_{1}_{2}'.format('roberts2017', netcdf, var))


def growingSeasonXR(region):
    """
    Input arrays should have the correct cropping already, so no need to trim in this function
    :param region: name of region
    :return:
    """

    prefix = 'roberts'
    crop = 'maize'          # Roberts only does maize, any other crop would require a new growing season definition

    if crop == 'winter_wheat':
        netcdf = 'wwh'
    elif crop == 'spring_wheat':
        netcdf = 'swh'
    else:
        netcdf = crop[:3]

    # plant = np.ma.zeros((dx, dy)) + 59  # March 1st
    # har = np.ma.zeros((dx, dy)) + 180   # 180 days growing season
    # ha = plant + har                    # add the 180 days to get the harvesting d.o.y.

    plant_date = '-03-01'
    har_date   = '-08-28'

    # The directory in which the climate data files are located
    f1 = '/project2/geos39650/ag_data/climate_projections/lpjml/'

    tasmin = xr.open_dataarray(f1 + 'tasmin_bced_1960_1999_hadgem2-es_rcp8p5_2005-2099_USA.nc4')
    tasmax = xr.open_dataarray(f1 + 'tasmax_bced_1960_1999_hadgem2-es_rcp8p5_2005-2099_USA.nc4')
    dx = tasmax.lat.size
    dy = tasmax.lon.size
    startyear = pd.to_datetime(tasmax.time.values[0]).year
    endyear   = pd.to_datetime(tasmax.time.values[-1]).year



    Gdd = np.ma.zeros((endyear - startyear + 1, dx, dy))
    Hdd = np.ma.zeros((endyear - startyear + 1, dx, dy))


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
        t_hrs = (t_hrs - 273.15)
        gdd_hrs = t_hrs.copy()

        t_high = 29
        t_low = 10
        ## GDD CALCULATION
        gdd_hrs[gdd_hrs > t_high] = t_high
        gdd_hrs = gdd_hrs - t_low
        gdd_hrs[gdd_hrs < 0] = 0
        Gdd[i, :, :] = np.sum(np.sum(gdd_hrs * (1 / 24), axis=3), axis=0)
        gdd_hrs = 0

        ## HDD CALCULATION
        hdd_hrs = t_hrs.copy() - t_high
        t_hrs = 0
        hdd_hrs[hdd_hrs < 0] = 0
        Hdd[i, :, :] = np.sum(np.sum(hdd_hrs * (1 / 24), axis=3), axis=0)
        hdd_hrs = 0

    Gdd[Gdd.mask] = 1e20
    growdir = '/project2/geos39650/ag_data/growseasons/'
    np.save(growdir + '{0}_{1}_gdd_{2}-{3}_{4}.npy'.format(prefix, netcdf, startyear, endyear, region), Gdd.data)
    print('gdd | saved')
    Hdd[Hdd.mask] = 1e20
    np.save(growdir + '{0}_{1}_hdd_{2}-{3}_{4}.npy'.format(prefix, netcdf, startyear, endyear, region), Hdd.data)
    print('hdd | saved')

    return True


def trueGS_OLD(temp_var, crop, model, t_shift, rainfed=True):
    """
    :param temp_var: CAN BE GDD, HDD, OR PR
    :param crop:
    :param model:
    :param t_shift:
    :param rainfed:
    :return:
    """
    netcdf = crop[:3]

    base_path = '/project2/ggcmi/AgMIP.output/{0}/phase2/{1}/A0'.format(model, crop)

    if temp_var == 'pr':
        temp_var_path = '/project2/geos39650/ag_data/agmerra/pr_agmerra_1980-2010.nc4'
        temp = xr.open_dataarray(temp_var_path).sel(lat=slice(49.25, 24.75), lon=slice(-124.75, -67.25))
    else:
        temp_var_path = '/project2/geos39650/ag_data/growseasons/{0}_USA_{1}.nc'.format(temp_var, t_shift)
        temp = xr.open_dataarray(temp_var_path)

    temp_lat     = temp.lat
    temp_lon     = temp.lon

    temp_out     = xr.DataArray(dims=['time', 'lat', 'lon'],
                               coords={'time':pd.date_range('1981', '2011', freq='1A'),
                                       'lat':temp_lat,
                                       'lon':temp_lon})

    w = 0 if rainfed else 'inf'
    N = 200 if model!='CARAIB' else 'NA'
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

        dayids = np.reshape(np.repeat(np.repeat(np.arange(1, lenyear2+1, 1), dx, axis=np.newaxis), dy, axis=np.newaxis),
                            (lenyear2, dx, dy))

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
                print("[{0}, {1}]".format(temp.lat.values[lat_i], temp.lon.values[lon_i]))
                pl.mask[0, lat_i, lon_i] = True     # Mask the earlier harvest, only show the latter
                ha.mask[0, lat_i, lon_i] = True     # Mask the earlier harvest, only show the latter

        # Collapse the 2 1-years into a single 2-year 2D array
        pl = np.sum(pl, axis=0)
        ha = np.sum(ha, axis=0)

        # Repeat out to fit daily size
        pl = np.repeat(pl[np.newaxis, :, :], lenyear2, axis=0)
        ha = np.repeat(ha[np.newaxis, :, :], lenyear2, axis=0)

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


def fixedGSlpjmlRCP(temp_var, historical):
    """
    :param temp_var: CAN BE GDD, HDD, OR PR
    :param crop:
    :param model:
    :param t_shift:
    :return:
    """
    crop = 'maize'
    netcdf = crop[:3]
    model = 'LPJmL'
    plant_date = '-03-01'
    har_date = '-08-28'

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
    else:
        temp_var_path = '/project2/geos39650/ag_data/growseasons/rcp_{0}_{1}_USA.nc'.format(temp_var, hist_ext)
        temp = xr.open_dataarray(temp_var_path)

    if historical:
        temp = temp.sel(time=slice('1951', '2004'))

    temp_lat     = temp.lat
    temp_lon     = temp.lon

    temp_out     = xr.DataArray(dims=['time', 'lat', 'lon'],
                               coords={'time':pd.date_range(str(startyear+1), str(endyear+1), freq='1A'),
                                       'lat':temp_lat,
                                       'lon':temp_lon})

    dy = endyear - startyear
    for i in range(1, dy + 1):
        year = startyear + i
        print(year)

        temp_i = temp.sel(time=slice(str(year) + plant_date, str(year) + har_date)).values
        temp_i = np.ma.masked_where(np.isnan(temp_i), temp_i)
        temp_out.loc[dict(time=str(year))] = np.ma.sum(temp_i, axis=0)

    save_dir = '/project2/geos39650/ag_data/growseasons/'
    temp_out.to_netcdf('{0}{1}_{2}_RCP_{3}_{4}_fixed.nc'.format(save_dir, model.lower(), netcdf, hist_ext, temp_var))
    return True


def trueGSlpjmlRCP_OLD(temp_var, historical, rainfed):
    """
    :param temp_var: CAN BE GDD, HDD, OR PR
    :param crop:
    :param model:
    :param t_shift:
    :param rainfed:
    :return:
    """
    crop = 'maize'
    netcdf = crop[:3]
    model = 'LPJmL'

    #C360_N200
    base_path = '/project2/ggcmi/AgMIP.output/LPJmL/phase2/fasttrack/HadGEM2-ES/rcp8p5/maize/c360_n200'

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
    else:
        temp_var_path = '/project2/geos39650/ag_data/growseasons/rcp_{0}_{1}_USA.nc'.format(temp_var, hist_ext)
        temp = xr.open_dataarray(temp_var_path)

    if historical:
        temp = temp.sel(time=slice('1951', '2004'))

    temp_lat     = temp.lat
    temp_lon     = temp.lon

    temp_out     = xr.DataArray(dims=['time', 'lat', 'lon'],
                               coords={'time':pd.date_range(str(startyear+1), str(endyear+1), freq='1A'),
                                       'lat':temp_lat,
                                       'lon':temp_lon})

    w = 'noirr2' if rainfed else 'firr2'
    plant_file  = "{0}/lpjml_HadGEM2-ES_rcp8p5_fullharm_plant-day_mai_global_annual_1951_2099_{1}.nc4".format(base_path, w)
    plant       = xr.open_dataarray(plant_file, decode_times=False).sel(lat=temp_lat, lon=temp_lon)
    plant['time'] = pd.date_range(str(1951), str(2100), freq='1A')
    plant       = plant.sel(time=slice(str(startyear), str(endyear)))
    plant       = np.ma.masked_array(plant)
    plant       = np.ma.masked_where(np.isnan(plant), plant)
    plant       = np.ma.harden_mask(plant).astype(int)

    har_file    = "{0}/lpjml_HadGEM2-ES_rcp8p5_fullharm_maty-day_mai_global_annual_1951_2099_{1}.nc4".format(base_path, w)
    har         = xr.open_dataarray(har_file, decode_times=False).sel(lat=temp_lat, lon=temp_lon)
    har['time'] = pd.date_range(str(1951), str(2100), freq='1A')
    har         = har.sel(time=slice(str(startyear), str(endyear)))
    har         = np.ma.masked_array(har)
    har         = np.ma.masked_where(np.isnan(har), har)
    har         = np.ma.harden_mask(har).astype(int)

    dy = endyear - startyear
    for i in range(0, dy):
        curryear = startyear + i
        print(curryear, '-', curryear + 1)

        temp_i      = temp.sel(time=slice(str(curryear), str(curryear + 1)))
        lenyear1    = temp_i.sel(time=str(curryear)).time.size
        lenyear2    = temp_i.time.size
        dx          = temp_i.lat.size
        dy          = temp_i.lon.size
        temp_i      = temp_i.values

        pl = plant[i:i + 2].copy()
        pl[1] = pl[1] + lenyear1
        ha = har[i:i + 2].copy()
        ha = pl + ha

        dayids = np.reshape(np.repeat(np.repeat(np.arange(0, lenyear2, 1), dx, axis=np.newaxis), dy, axis=np.newaxis),
                            (lenyear2, dx, dy))

        # Cut harvest dates before the year started
        pl = np.ma.masked_where(ha <= lenyear1, pl)
        ha = np.ma.masked_where(ha <= lenyear1, ha)
        # Cut harvest dates after the year ended
        pl = np.ma.masked_where(ha >= lenyear2, pl)
        ha = np.ma.masked_where(ha >= lenyear2, ha)

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
                print("[{0}, {1}]".format(temp.lat.values[lat_i], temp.lon.values[lon_i]))
                pl.mask[0, lat_i, lon_i] = True     # Mask the latter harvest, only show the previous
                ha.mask[0, lat_i, lon_i] = True     # Mask the latter harvest, only show the previous


        # Collapse the 2 1-years into a single 2-year 2D array
        pl = np.sum(pl, axis=0)
        ha = np.sum(ha, axis=0)

        # Repeat out to fit daily size
        pl = np.repeat(pl[np.newaxis, :, :], lenyear2, axis=0)
        ha = np.repeat(ha[np.newaxis, :, :], lenyear2, axis=0)

        # Assert that same grid cells are masked in PL and HA arrays
        assert np.array_equal(ha.mask, pl.mask), "Masks of dates don't match for PL and HA"

        temp_i = np.ma.masked_where(dayids < pl, temp_i)
        temp_i = np.ma.masked_where(dayids > ha, temp_i)
        temp_i = np.ma.masked_where(temp_i >= 1e20, temp_i)
        temp_i = np.ma.harden_mask(temp_i)

        temp_out.loc[dict(time=str(curryear + 1))] = np.ma.sum(temp_i, axis=0)

    save_dir = '/project2/geos39650/ag_data/true_gs/'
    temp_out.to_netcdf('{0}{1}_{2}_RCP_{3}_{4}_{5}_test.nc'.format(save_dir, model.lower(), netcdf, hist_ext, temp_var, w))
    return True


def trueGS(temp_var, crop, model, t_shift, rainfed=True):
    """
    :param temp_var: CAN BE GDD, HDD, OR PR
    :param crop:
    :param model:
    :param t_shift:
    :param rainfed:
    :return:
    """
    netcdf = crop[:3]
    startyear = 1981
    endyear   = 2011

    base_path = '/project2/ggcmi/AgMIP.output/{0}/phase2/{1}/A0'.format(model, crop)

    if temp_var == 'pr':
        temp_var_path = '/project2/geos39650/ag_data/agmerra/pr_agmerra_1980-2010.nc4'
        temp = xr.open_dataarray(temp_var_path).sel(lat=slice(49.25, 24.75), lon=slice(-124.75, -67.25))
    else:
        temp_var_path = '/project2/geos39650/ag_data/growseasons/{0}_USA_{1}.nc'.format(temp_var, t_shift)
        temp = xr.open_dataarray(temp_var_path)

    temp = temp.sel(time=slice('1981', '2010'))

    temp_lat = temp.lat
    temp_lon = temp.lon

    extra_year = xr.DataArray(dims=['time', 'lat', 'lon'],
                            coords={'time': pd.date_range(str(endyear)+'-01-01', str(endyear)+'-12-31', freq='1D'),
                                    'lat': temp_lat,
                                    'lon': temp_lon})

    temp = xr.concat([temp, extra_year], dim='time')

    temp_out     = xr.DataArray(dims=['time', 'lat', 'lon'],
                               coords={'time':pd.date_range(str(startyear), str(endyear), freq='1A'),
                                       'lat':temp_lat,
                                       'lon':temp_lon})

    w = 0 if rainfed else 'inf'
    N = 200 if model!='CARAIB' else 'NA'
    plant_file  = "{0}/plant-day/{1}_agmerra_fullharm_plant-day_mai_global_annual_1980_2010_C360_T{2}_W{3}_N{4}_A0.nc4".format(base_path, model.lower(), t_shift, w, N)
    plant       = xr.open_dataarray(plant_file, decode_times=False).sel(lat=temp_lat, lon=temp_lon)
    plant['time'] = pd.date_range(str(1980), str(2011), freq='1A')
    plant       = plant.sel(time=slice(str(startyear), str(endyear-1))).values
    plant       = np.append(plant, np.empty(plant[0].shape)[np.newaxis,:, :]*np.NaN, axis=0)
    plant       = np.ma.masked_array(plant)
    plant       = np.ma.masked_where(np.isnan(plant), plant)
    plant       = np.ma.harden_mask(plant).astype(int)

    har_file    = "{0}/maty-day/{1}_agmerra_fullharm_maty-day_mai_global_annual_1980_2010_C360_T{2}_W{3}_N{4}_A0.nc4".format(base_path, model.lower(), t_shift, w, N)
    har         = xr.open_dataarray(har_file, decode_times=False).sel(lat=temp_lat, lon=temp_lon)
    har['time'] = pd.date_range(str(1980), str(2011), freq='1A')
    har         = har.sel(time=slice(str(startyear), str(endyear-1))).values
    har         = np.append(har, np.empty(har[0].shape)[np.newaxis, :, :]*np.NaN, axis=0)
    har         = np.ma.masked_array(har)
    har         = np.ma.masked_where(np.isnan(har), har)
    har         = np.ma.harden_mask(har).astype(int)

    dy = endyear - startyear
    for i in range(0, dy):
        curryear = startyear + i
        print(curryear, '-', curryear + 1)

        temp_i = temp.sel(time=slice(str(curryear), str(curryear + 1)))
        lenyear1 = temp_i.sel(time=str(curryear)).time.size
        lenyear2 = temp_i.time.size
        dx = temp_i.lat.size
        dy = temp_i.lon.size
        temp_i = temp_i.values

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

        dayids = np.reshape(
            np.repeat(np.repeat(np.arange(1, lenyear2 + 1, 1), dx, axis=np.newaxis), dy, axis=np.newaxis),
            (lenyear2, dx, dy))

        # Mask PLANTING dates AFTER Year-1 ended
        pl = np.ma.masked_where(pl > lenyear1, pl)
        ha = np.ma.masked_where(pl > lenyear1, ha)

        # Mask HARVEST dates AFTER Year-2 ended
        pl = np.ma.masked_where(ha > lenyear2, pl)
        ha = np.ma.masked_where(ha > lenyear2, ha)

        if i == dy - 1:  # If it's the last year of calculation
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
                print("[{0}, {1}]".format(temp.lat.values[lat_i], temp.lon.values[lon_i]))
                pl.mask[1, lat_i, lon_i] = True  # Mask the earlier harvest, only show the latter
                ha.mask[1, lat_i, lon_i] = True  # Mask the earlier harvest, only show the latter

        # Collapse the 2 1-years into a single 2-year 2D array
        pl = np.sum(pl, axis=0)
        ha = np.sum(ha, axis=0)

        # Repeat out to fit daily size
        pl = np.repeat(pl[np.newaxis, :, :], lenyear2, axis=0)
        ha = np.repeat(ha[np.newaxis, :, :], lenyear2, axis=0)

        # Assert that same grid cells are masked in PL and HA arrays
        assert np.array_equal(ha.mask, pl.mask), "Masks of dates don't match for PL and HA"

        temp_i = np.ma.masked_where(dayids < pl, temp_i)
        temp_i = np.ma.masked_where(dayids > ha, temp_i)
        temp_i = np.ma.masked_where(temp_i >= 1e20, temp_i)
        temp_i = np.ma.harden_mask(temp_i)

        temp_out.loc[dict(time=str(curryear))] = np.ma.sum(temp_i, axis=0)

    save_dir = '/project2/geos39650/ag_data/true_gs/'
    temp_out.to_netcdf('{0}{1}_{2}_TrueGS_{3}_T{4}_W{5}.nc'.format(save_dir, model.lower(), netcdf, temp_var, t_shift, w))
    return True


def trueGSlpjmlRCP(temp_var, historical, rainfed):
    """
    :param temp_var: CAN BE GDD, HDD, OR PR
    :param crop:
    :param model:
    :param t_shift:
    :param rainfed:
    :return:
    """
    crop = 'maize'
    netcdf = crop[:3]
    model = 'LPJmL'

    #C360_N200
    base_path = '/project2/ggcmi/AgMIP.output/LPJmL/phase2/fasttrack/HadGEM2-ES/rcp8p5/maize/c360_n200'

    if historical:
        hist_ext = 'historical_1950-2004'
        startyear   = 1952 # b/c the pl/har files start with 1951
        endyear     = 2005
        save_ext = 'historical_1952-2004'
    else:
        hist_ext = 'rcp8p5_2005-2099'
        startyear   = 2006
        endyear     = 2100
        save_ext = 'rcp8p5_2006-2099'

    if temp_var == 'pr':
        temp_var_path = '/project2/geos39650/ag_data/climate_projections/lpjml/pr_bced_1960_1999_hadgem2-es_{0}_USA.nc4'.format(hist_ext)
        temp = xr.open_dataarray(temp_var_path)
    else:
        temp_var_path = '/project2/geos39650/ag_data/growseasons/rcp_{0}_{1}_USA.nc'.format(temp_var, hist_ext)
        temp = xr.open_dataarray(temp_var_path)

    if historical:
        temp = temp.sel(time=slice('1952', '2004'))
    else:
        temp = temp.sel(time=slice('2006', '2099'))

    temp_lat = temp.lat
    temp_lon = temp.lon

    extra_year = xr.DataArray(dims=['time', 'lat', 'lon'],
                            coords={'time': pd.date_range(str(endyear)+'-01-01', str(endyear)+'-12-31', freq='1D'),
                                    'lat': temp_lat,
                                    'lon': temp_lon})

    temp = xr.concat([temp, extra_year], dim='time')

    temp_out     = xr.DataArray(dims=['time', 'lat', 'lon'],
                                coords={'time':pd.date_range(str(startyear), str(endyear), freq='1A'),
                                       'lat':temp_lat,
                                       'lon':temp_lon})

    w = 'noirr2' if rainfed else 'firr2'
    plant_file  = "{0}/lpjml_HadGEM2-ES_rcp8p5_fullharm_plant-day_mai_global_annual_1951_2099_{1}.nc4".format(base_path, w)
    plant       = xr.open_dataarray(plant_file, decode_times=False).sel(lat=temp_lat, lon=temp_lon)
    plant['time'] = pd.date_range(str(1951), str(2100), freq='1A')
    plant       = plant.sel(time=slice(str(startyear), str(endyear-1))).values
    plant       = np.append(plant, np.empty(plant[0].shape)[np.newaxis,:, :]*np.NaN, axis=0)
    plant       = np.ma.masked_array(plant)
    plant       = np.ma.masked_where(np.isnan(plant), plant)
    plant       = np.ma.harden_mask(plant).astype(int)

    har_file    = "{0}/lpjml_HadGEM2-ES_rcp8p5_fullharm_maty-day_mai_global_annual_1951_2099_{1}.nc4".format(base_path, w)
    har         = xr.open_dataarray(har_file, decode_times=False).sel(lat=temp_lat, lon=temp_lon)
    har['time'] = pd.date_range(str(1951), str(2100), freq='1A')
    har         = har.sel(time=slice(str(startyear), str(endyear-1))).values
    har         = np.append(har, np.empty(har[0].shape)[np.newaxis, :, :]*np.NaN, axis=0)
    har         = np.ma.masked_array(har)
    har         = np.ma.masked_where(np.isnan(har), har)
    har         = np.ma.harden_mask(har).astype(int)

    dy = endyear - startyear
    for i in range(0, dy):
        curryear = startyear + i
        print(curryear, '-', curryear + 1)

        temp_i      = temp.sel(time=slice(str(curryear), str(curryear + 1)))
        lenyear1    = temp_i.sel(time=str(curryear)).time.size
        lenyear2    = temp_i.time.size
        dx          = temp_i.lat.size
        dy          = temp_i.lon.size
        temp_i      = temp_i.values

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
                print("[{0}, {1}]".format(temp.lat.values[lat_i], temp.lon.values[lon_i]))
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

        temp_i = np.ma.masked_where(dayids < pl, temp_i)
        temp_i = np.ma.masked_where(dayids > ha, temp_i)
        temp_i = np.ma.masked_where(temp_i >= 1e20, temp_i)
        temp_i = np.ma.harden_mask(temp_i)

        temp_out.loc[dict(time=str(curryear))] = np.ma.sum(temp_i, axis=0)

    save_dir = '/project2/geos39650/ag_data/true_gs/'
    temp_out.to_netcdf('{0}{1}_{2}_RCP_{3}_{4}_{5}.nc'.format(save_dir, model.lower(), netcdf, save_ext, temp_var, w))
    return True


def RCPgs():
    crop = 'maize'
    netcdf = crop[:3]

    # get GDD and HDD
    growingSeasonXR('USA')


    # # get pr values
    # var = 'pr'
    # aggfxn = lambda x: np.ma.sum(x, axis=0)
    # save_name = '{0}_{1}_pr_{2}-{3}_{4}'.format('roberts', netcdf, 2031, 2099, 'USA')
    # growingSeasonVarXR(var, aggfxn, '-03-01', '-08-28', save_name)


if __name__ == '__main__':
    temp_vars = ['gdd', 'hdd', 'pr']
    # for var in temp_vars:
    #     his = trueGSlpjmlRCP(var, historical=True, rainfed=True)
    #     fut = trueGSlpjmlRCP(var, historical=False, rainfed=True)

    crop     = 'maize'
    models   = ['LPJmL', 'LPJ-GUESS', 'CARAIB', 'PEPIC', 'GEPIC', 'EPIC-TAMU', 'pDSSAT']
    t_shifts = np.arange(-1,7)

    for var in temp_vars:
        for model in models:
            for t_shift in t_shifts:
                try:
                    outfile = trueGS(var, crop, model, t_shift, rainfed=True)
                except:
                    print('NOT FOUND: {0}, {1}, T+{2}'.format(model, var, t_shift))




