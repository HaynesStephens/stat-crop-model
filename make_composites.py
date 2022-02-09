import xarray as xr
import numpy as np

def loadGDDandHDD():
    load_path = '/project2/geos39650/ag_data/agmerra/'
    save_path = '/project2/geos39650/ag_data/growseasons/'
    loc = [81, 131, 110, 226]
    top, bottom, left, right = loc
    years = np.arange(1980,2011)


    tasmin = xr.open_dataarray(load_path + 'tasmin' + '_agmerra_1980-2010.nc4').isel(lat=slice(top,bottom),
                                                                               lon=slice(left,right))
    tasmax = xr.open_dataarray(load_path + 'tasmax' + '_agmerra_1980-2010.nc4').isel(lat=slice(top,bottom),
                                                                               lon=slice(left,right))

    for t_shift in [-1, 0, 1, 2, 3, 4, 5, 6]:
        gdd = xr.DataArray(np.zeros(tasmin.shape),
                           coords={'time': tasmin.time, 'lat': tasmin.lat, 'lon': tasmin.lon},
                           dims=['time', 'lat', 'lon'])
        hdd = xr.DataArray(np.zeros(tasmin.shape),
                           coords={'time': tasmin.time, 'lat': tasmin.lat, 'lon': tasmin.lon},
                           dims=['time', 'lat', 'lon'])

        for y in years:
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
            t_hrs = t_hrs - 273.15 + t_shift
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

        gdd.to_netcdf('{0}gdd_USA_{1}.nc'.format(save_path, t_shift))
        hdd.to_netcdf('{0}hdd_USA_{1}.nc'.format(save_path, t_shift))

    return True


def loadGDDandHDD_lpjml(historical=False):
    load_path = '/project2/geos39650/ag_data/climate_projections/lpjml/'
    save_path = '/project2/geos39650/ag_data/growseasons/'
    # loc = [81, 131, 110, 226]
    # top, bottom, left, right = loc
    if historical:
        years = np.arange(1950,2005)
        file_ext = '_bced_1960_1999_hadgem2-es_historical_1950-2004_USA.nc4'
        save_ext = 'historical_1950-2004_USA.nc4'
    else:
        years = np.arange(2005,2100)
        file_ext = '_bced_1960_1999_hadgem2-es_rcp8p5_2005-2099_USA.nc4'
        save_ext = 'rcp8p5_2005-2099_USA.nc4'


    tasmin = xr.open_dataarray(load_path + 'tasmin' + file_ext)
    tasmax = xr.open_dataarray(load_path + 'tasmax' + file_ext)

    gdd = xr.DataArray(np.zeros(tasmin.shape),
                       coords={'time': tasmin.time, 'lat': tasmin.lat, 'lon': tasmin.lon},
                       dims=['time', 'lat', 'lon'])
    hdd = xr.DataArray(np.zeros(tasmin.shape),
                       coords={'time': tasmin.time, 'lat': tasmin.lat, 'lon': tasmin.lon},
                       dims=['time', 'lat', 'lon'])

    for y in years:
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

    gdd.to_netcdf('{0}rcp_gdd_{1}'.format(save_path, save_ext))
    hdd.to_netcdf('{0}rcp_hdd_{1}'.format(save_path, save_ext))

    return True


def loadBins():
    load_path = '/project2/geos39650/ag_data/agmerra/'
    save_path = '/project2/geos39650/ag_data/growseasons/'
    loc = [81, 131, 110, 226]
    top, bottom, left, right = loc
    years = np.arange(1980,2011)


    tasmin = xr.open_dataarray(load_path + 'tasmin' + '_agmerra_1980-2010.nc4').isel(lat=slice(top,bottom),
                                                                               lon=slice(left,right))
    tasmax = xr.open_dataarray(load_path + 'tasmax' + '_agmerra_1980-2010.nc4').isel(lat=slice(top,bottom),
                                                                               lon=slice(left,right))

    # for t_shift in [-1, 0, 1, 2, 3, 4, 5, 6]:
    for t_shift in [6]:
        tbins = xr.DataArray(np.zeros((tasmin.shape[0], tasmin.shape[1], tasmin.shape[2], 14)),
                           coords={'time': tasmin.time, 'lat': tasmin.lat, 'lon': tasmin.lon, 'tbin': np.arange(14)},
                           dims=['time', 'lat', 'lon', 'tbin'])

        for y in years:
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
            t_hrs = t_hrs - 273.15 + t_shift

            t_hrs = np.clip(t_hrs, a_min=t_hrs.min(), a_max=42)

            bins = np.arange(0, 43, 3)
            temp = np.apply_along_axis(lambda a: np.histogram(a, bins=bins)[0], 3, t_hrs) / 24

            tbins.loc[dict(time=tasmini.time)] = temp

        tbins.to_netcdf('{0}tbins_USA_{1}.nc'.format(save_path, t_shift))

    return True


def loadBins_lpjml(historical=False):
    load_path = '/project2/geos39650/ag_data/climate_projections/lpjml/'
    save_path = '/project2/geos39650/ag_data/growseasons/'
    # loc = [81, 131, 110, 226]
    # top, bottom, left, right = loc
    if historical:
        years = np.arange(1950,2005)
        file_ext = '_bced_1960_1999_hadgem2-es_historical_1950-2004_USA.nc4'
        save_ext = 'historical_1950-2004_USA.nc4'
    else:
        years = np.arange(2005,2100)
        file_ext = '_bced_1960_1999_hadgem2-es_rcp8p5_2005-2099_USA.nc4'
        save_ext = 'rcp8p5_2005-2099_USA.nc4'


    tasmin = xr.open_dataarray(load_path + 'tasmin' + file_ext)
    tasmax = xr.open_dataarray(load_path + 'tasmax' + file_ext)

    tbins = xr.DataArray(np.zeros((tasmin.shape[0], tasmin.shape[1], tasmin.shape[2], 14)),
                         coords={'time': tasmin.time, 'lat': tasmin.lat, 'lon': tasmin.lon, 'tbin': np.arange(14)},
                         dims=['time', 'lat', 'lon', 'tbin'])

    for y in years:
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

        t_hrs = np.clip(t_hrs, a_min=t_hrs.min(), a_max=42)

        bins = np.arange(0, 43, 3)
        temp = np.apply_along_axis(lambda a: np.histogram(a, bins=bins)[0], 3, t_hrs) / 24

        tbins.loc[dict(time=tasmini.time)] = temp

    tbins.to_netcdf('{0}rcp_tbins_{1}.nc'.format(save_path, save_ext))

    return True



if __name__ == '__main__':
    first  = loadBins_lpjml(historical=True)
    second = loadBins_lpjml(historical=False)
