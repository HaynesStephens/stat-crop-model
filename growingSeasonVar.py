import xarray as xr
import numpy as np
import pandas as pd

def growingSeasonVar(var, aggfxn, pl_in, ha_in, save_name, loc=None):
    """
    A function that can calculate aggregated growing-season values from the AgMERRA dataset
    :param var: The variable you want (i.e. 'tas', 'tasmin', 'tasmax', 'pr')
    :param aggfxn: The aggregating function. How you want to aggregate to a single growing-season value.
           Note: This function will be performed on a cut-out array from the AgMERRA data.
                 That cut-out array will have dimensions of (730, 360, 720), a 2-year span, in order to account
                 for growing seasons that span between two calendar years.
    :param pl: An array of your planting days, in units of "Day of Year". This needs to be a 2D array of dimension (360, 720).
           Note: This means that planting days for a given location will be the same for every growing season.
                 This applies to harvest day as well.
    :param har: Same as planting day array, except for harvesting days.
    :param save_name: The name of the file you want to save, excluding the file extension. Files will be saved as '.npy' arrays.
    :param loc: A list specifying the specific x,y area you want to cover.
                LOC is a 4-item list of the form: [xStart, xEnd, yStart, yEnd], so you'll need to know the indices of the
                region you want. Also, remember that Python will exclude the final index value in the cut-out array.
    :return: Save a '.npy' file of shape (30, 360, 720) with the appropriate growing-season values from the year 1981 to 2010.
             Note: growing-season values for 1981 will correspond to the yields harvested in 1981, and so on.
    """
    # If LOC is specified, cut down the AgMERRA data to the appropriate area.
    if loc != None:
        top, bottom, left, right = loc
    else:
        top, bottom, left, right = [0, 360, 0, 720]
    dx = bottom - top
    dy = right - left
    # Make copies of the plant and harvest arrays, for safe measure.
    # Broadcast the 2D array into a 3D array of length 730 days. This serves as a 2-year window of plant/harvest days.
    pl  = np.repeat(pl_in.copy()[np.newaxis, :, :], 730, axis=0)
    har = np.repeat(ha_in.copy()[np.newaxis, :, :], 730, axis=0)
    pl  = pl[:, top:bottom, left:right]
    har = har[:, top:bottom, left:right]
    # Shift planting and harvest days so that all harvest takes place in the final year of the 2-year window.
    pl  = np.where(har < 365, pl + 365, pl)
    har = np.where(har < 365, har + 365, har)
    # Create an array that gives day-of-year, to be used in order to cut-out growing seasons.
    dayids = np.reshape(np.repeat(np.repeat(np.arange(0, 730, 1), dx, axis=np.newaxis), dy, axis=np.newaxis), (730, dx, dy))

    # The directory in which the AgMERRA data files are located
    f1 = '/project2/geos39650/ag_data/agmerra/'

    # Create an empty 30-year array that will be filled with calculated growing-season values.
    var_out = np.ma.zeros((30, dx, dy))

    # The first year of growing-season values to be calculated, this corresponds to i=0 in the following loop.
    startyear = 1981
    # The first leap year encountered in the loop, following leaps are +4
    leapyear = 1984
    # Add a day for each leap year, cumulatively through the time loop
    leap = 0

    # Loop through the AgMERRA set using a 2-year window and determine growing-season values.
    for i in range(30):
        print(startyear + i)
        # If a leap year is encountered, add a day to the indexing variable.
        if (i + 1) % 4 == 0:
            leap += 1
            print('leap year')

        # Load the relevant variable from the AgMERRA dataset. Cut out the 2-year window that corresponds to a particular growing season
        var_arri    = hfile(f1 + var + '_agmerra_1980-2010.nc4', 'r')[var][((i) * 365 + leap):((i) * 365 + 730 + leap), top:bottom, left:right]

        # Mask that 2-year window outside of the growing seaons,
        # i.e. where day-of-year is less than plant day and greater than harvest day.
        var_arri = np.ma.masked_where(dayids < pl, var_arri)
        var_arri = np.ma.masked_where(dayids > har, var_arri)
        # Mask the array where values are fill-in, i.e. n/a values, denoted in this dataset as 1e20
        var_arri = np.ma.masked_where(var_arri >= 1e20, var_arri)
        # Harden the masked array so things don't act funny (they still might, but oh well).
        var_arri = np.ma.harden_mask(var_arri)
        # Apply the AGGFXN to the daily values in the growing season, and add these values into the outputted array
        # at the appropriate year index.
        var_out[i, :, :] = aggfxn(var_arri)

    # Remask the final outputted array with 1e20 values and save as a .npy file under SAVE_NAME
    var_out[var_out.mask] = 1e20
    growdir = '/project2/geos39650/ag_data/growseasons/'
    np.save(growdir + '{0}.npy'.format(save_name), var_out.data)
    print('{0} | saved'.format(save_name))
    return var_out

"""EXAMPLE RUN
growingSeasonVar('tas', lambda x: np.mean(x, axis=0), np.ones((360, 720))*100, np.ones((360, 720))*280, 'test')
"""

def growingSeasonVarXR(var, aggfxn, plant_date, har_date, save_name):
    """
    A function that can calculate aggregated growing-season values from the AgMERRA dataset
    :param var: The variable you want (i.e. 'tas', 'tasmin', 'tasmax', 'pr')
    :param aggfxn: The aggregating function. How you want to aggregate to a single growing-season value.
           Note: This function will be performed on a cut-out array from the AgMERRA data.
                 That cut-out array will have dimensions of (730, 360, 720), a 2-year span, in order to account
                 for growing seasons that span between two calendar years.
    :param pl: An array of your planting days, in units of "Day of Year". This needs to be a 2D array of dimension (360, 720).
           Note: This means that planting days for a given location will be the same for every growing season.
                 This applies to harvest day as well.
    :param har: Same as planting day array, except for harvesting days.
    :param save_name: The name of the file you want to save, excluding the file extension. Files will be saved as '.npy' arrays.
    :param loc: A list specifying the specific x,y area you want to cover.
                LOC is a 4-item list of the form: [xStart, xEnd, yStart, yEnd], so you'll need to know the indices of the
                region you want. Also, remember that Python will exclude the final index value in the cut-out array.
    :return: Save a '.npy' file of shape (30, 360, 720) with the appropriate growing-season values from the year 1981 to 2010.
             Note: growing-season values for 1981 will correspond to the yields harvested in 1981, and so on.
    """

    # The directory in which the climate data files are located
    f1 = '/project2/geos39650/ag_data/climate_projections/lpjml/'

    var_arr = xr.open_dataarray(f1 + var + '_bced_1960_1999_hadgem2-es_rcp8p5_2031-2099_USA.nc4')
    dx = var_arr.lat.size
    dy = var_arr.lon.size
    startyear = pd.to_datetime(var_arr.time.values[0]).year
    endyear   = pd.to_datetime(var_arr.time.values[-1]).year

    # Create an empty 30-year array that will be filled with calculated growing-season values.
    var_out = np.ma.zeros((endyear - startyear + 1, dx, dy))

    # Loop through the AgMERRA set using a 2-year window and determine growing-season values.
    for i in range(endyear - startyear + 1):
        year = startyear + i
        print(year)

        # Load the relevant variable from the AgMERRA dataset. Cut out the time window that corresponds to a particular growing season
        var_arri = var_arr.sel(time=slice(str(year) + plant_date, str(year) + har_date)).values
        var_arri = np.ma.masked_where(np.isnan(var_arri), var_arri)

        # Harden the masked array so things don't act funny (they still might, but oh well).
        var_arri = np.ma.harden_mask(var_arri)

        # Apply the AGGFXN to the daily values in the growing season, and add these values into the outputted array
        # at the appropriate year index.
        var_out[i, :, :] = aggfxn(var_arri)

    # Remask the final outputted array with 1e20 values and save as a .npy file under SAVE_NAME
    var_out[var_out.mask] = 1e20
    growdir = '/project2/geos39650/ag_data/growseasons/'
    np.save(growdir + '{0}.npy'.format(save_name), var_out.data)
    print('{0} | saved'.format(save_name))
    return var_out

"""EXAMPLE RUN
growingSeasonVarXR('pr', lambda x: np.sum(x, axis=0), '-03-01', '-08-28', 'test')
"""


# def growingSeasonQuant(vars, aggfxn, pl_in, ha_in, save_name, loc=None):
#     """
#     A function that can calculate aggregated growing-season values from the AgMERRA dataset
#     :param vars: The LIST of variables you want (i.e. ['tas', 'tasmin', 'tasmax', 'pr'])
#     :param aggfxn: The aggregating function that works on the LIST of variable arrays to
#                     get your growing season values, i.e. How you want to aggregate to a single growing-season value.
#            Note: This function will be performed on a cut-out array from the AgMERRA data.
#                  That cut-out array will have dimensions of (730, 360, 720), a 2-year span, in order to account
#                  for growing seasons that span between two calendar years.
#     :param pl: A 3D array of your planting days, in units of "Day of Year". This needs to be a 3D array of dimension (30, 360, 720).
#            Note: Should only include years for yields, i.e. 1981-2010.
#            This allows for year-to-year variation in plant/harvest days.
#     :param har: planting day PLUS MATY DAY.
#     :param save_name: The name of the file you want to save, excluding the file extension. Files will be saved as '.npy' arrays.
#     :param loc: A list specifying the specific x,y area you want to cover.
#                 LOC is a 4-item list of the form: [xStart, xEnd, yStart, yEnd], so you'll need to know the indices of the
#                 region you want. Also, remember that Python will exclude the final index value in the cut-out array.
#     :return: Save a '.npy' file of shape (30, 360, 720) with the appropriate growing-season values from the year 1981 to 2010.
#              Note: growing-season values for 1981 will correspond to the yields harvested in 1981, and so on.
#     """
#     # If LOC is specified, cut down the AgMERRA data to the appropriate area.
#     if loc != None:
#         top, bottom, left, right = loc
#     else:
#         top, bottom, left, right = [0, 360, 0, 720]
#     dx = bottom - top
#     dy = right - left
#
#     # Create an array that gives day-of-year, to be used in order to cut-out growing seasons.
#     dayids = np.reshape(np.repeat(np.repeat(np.arange(0, 730, 1), dx, axis=np.newaxis), dy, axis=np.newaxis), (730, dx, dy))
#
#     # The directory in which the AgMERRA data files are located
#     f1 = '/project2/geos39650/ag_data/agmerra/'
#
#     var_out = np.ma.zeros((30, dx, dy))
#
#     # The first year of growing-season values to be calculated, this corresponds to i=0 in the following loop.
#     startyear = 1981
#     # The first leap year encountered in the loop, following leaps are +4
#     leapyear = 1984
#     # Add a day for each leap year, cumulatively through the time loop
#     leap = 0
#
#     # Loop through the AgMERRA set using a 2-year window and determine growing-season values.
#     for i in range(30):
#         print(startyear + i)
#         # If a leap year is encountered, add a day to the indexing variable.
#         if (i + 1) % 4 == 0:
#             leap += 1
#             print('leap year')
#
#         # Make copies of the plant and harvest arrays, for safe measure.
#         # Broadcast the 2D array into a 3D array of length 730 days. This serves as a 2-year window of plant/harvest days.
#         pl = pl_in[i:i+2, top:bottom, left:right]
#         har = ha_in[i:i+2, top:bottom, left:right]
#         # Broadcast the 2D array into a 3D array of length 730 days. This serves as a 2-year window of plant/harvest days.
#         pl = np.repeat(pl[np.newaxis, :, :], 730, axis=0)
#         har = np.repeat(har[np.newaxis, :, :], 730, axis=0)
#         # Shift planting and harvest days so that all harvest takes place in the final year of the 2-year window.
#         pl = np.where(har < 365, pl + 365, pl)
#         har = np.where(har < 365, har + 365, har)
#
#         vars_in = []
#         for k in range(len(vars)):
#             var = vars[k]
#             # Load the relevant variable from the AgMERRA dataset. Cut out the 2-year window that corresponds to a particular growing season
#             var_arri    = hfile(f1 + var + '_agmerra_1980-2010.nc4', 'r')[var][((i) * 365 + leap):((i) * 365 + 730 + leap), top:bottom, left:right]
#
#             # Mask that 2-year window outside of the growing seaons,
#             # i.e. where day-of-year is less than plant day and greater than harvest day.
#             var_arri = np.ma.masked_where(dayids < pl, var_arri)
#             var_arri = np.ma.masked_where(dayids > har, var_arri)
#             # Mask the array where values are fill-in, i.e. n/a values, denoted in this dataset as 1e20
#             var_arri = np.ma.masked_where(var_arri >= 1e20, var_arri)
#             # Harden the masked array so things don't act funny (they still might, but oh well).
#             var_arri = np.ma.harden_mask(var_arri)
#             # Apply the AGGFXN to the daily values in the growing season, and add these values into the outputted array
#             # at the appropriate year index.
#             vars_in[k] = var_arri
#             vars_in[k] = np.ma.harden_mask(vars_in[k])
#         var_out[i, :, :] = aggfxn(vars_in)
#
#     # Remask the final outputted array with 1e20 values and save as a .npy file under SAVE_NAME
#     var_out[var_out.mask] = 1e20
#     growdir = '/project2/geos39650/ag_data/growseasons/'
#     np.save(growdir + '{0}.npy'.format(save_name), var_out.data)
#     print('{0} | saved'.format(save_name))
#     return var_out
#
# """EXAMPLE RUN
# growingSeasonQuant(['tasmin', 'tasmax'], lambda lst: (np.mean(lst[0], axis=0) + np.mean(lst[1], axis=0))/2,
#                     np.ones((730, 360, 720))*100, np.ones((730, 360, 720))*280, 'quant_test', loc = [82, 108, 152, 200])
# """
