#!/bin/env python
from h5py import File as hfile
import numpy as np
from cropArea import getarea


def loadData(name):
    f1 = '/project2/geos39650/ag_data/agmerra/'
    print(name + ' | loading')
    arr = hfile(f1 + name + '_agmerra_1980-2010.nc4', 'r')[name][:]
    print(name + ' | loaded')
    return arr


#######################################################################################################################

def growingSeasonAll(tas, tasmax, tasmin, prec, season_type = 'model', model = 'LPJmL'):
    print('START')
    ### Aggregate a gridded climate var to growing season ###
    # * args * #
    # var: input array
    dayids = np.reshape(np.repeat(np.repeat(np.arange(0, 730, 1), 360, axis=np.newaxis), 720, axis=np.newaxis),
                        (730, 360, 720))

    # print('tas | loading')
    # tas = hfile(f1 + 'tas_agmerra_1980-2010.nc4', 'r')['tas'][:]
    # print('tas | loaded')
    #
    # print('pr | loading')
    # prec = hfile(f1 + 'pr_agmerra_1980-2010.nc4', 'r')['pr'][:]
    # print('pr | loaded')

    crops = ['maize', 'rice', 'soy', 'winter_wheat', 'spring_wheat']

    def clipAndMask(arr_input, pl, ha, dayids, leap):
        arr_i = arr_input[((i) * 365 + leap):((i) * 365 + 730 + leap), :, :]
        arr_i = np.ma.masked_where(dayids < pl, arr_i)
        arr_i = np.ma.masked_where(dayids > ha, arr_i)
        return arr_i

    for crop in crops:
        print('crop | %s' % crop)
        if season_type == 'model':
            pl, ha = modelGrowingSeason(model, crop)
        elif season_type == 'LobellField':
            pl, ha = LobellFieldGrowingSeason(crop)
            rm, rmI, CAL, netcdf, nvar = getarea(crop)
            area = rm + rmI
            area[area==1e20] = 0

        tas_out     = np.zeros((30, 360, 720))
        tasmin_out  = np.zeros((30, 360, 720))
        tasmax_out  = np.zeros((30, 360, 720))
        pr_out      = np.zeros((30, 360, 720))

        startyear = 1981 # this is index 0 in the following loop
        leapyear = 1984 # this is the first leap year, following leaps are +4
        leap = 0 # add a day for each leap year cumulatively through the time loop
        for i in range(30):
            print(startyear + i)
            if (i + 1) % 4 == 0:
                leap += 1
                print('leap year')
            # tasi = tas[((i) * 365 + leap):((i) * 365 + 730 + leap), :, :]
            # tasi = np.ma.masked_where(dayids < pl, tasi)
            # tasi = np.ma.masked_where(dayids > ha, tasi)
            # pri = prec[((i) * 365 + leap):((i) * 365 + 730 + leap), :, :]
            # pri = np.ma.masked_where(dayids < pl, pri)
            # pri = np.ma.masked_where(dayids > ha, pri)

            tasi    = clipAndMask(tas, pl, ha, dayids, leap)
            tasmini = clipAndMask(tasmin, pl, ha, dayids, leap)
            tasmaxi = clipAndMask(tasmax, pl, ha, dayids, leap)
            pri     = clipAndMask(prec, pl, ha, dayids, leap)

            #Gathers mean growing-season values, can change for different growing season values
            if season_type == 'model':
                tas_out[i, :, :]    = np.ma.mean(tasi, axis=0)
                tasmin_out[i, :, :] = np.ma.mean(tasmini, axis=0)
                tasmax_out[i, :, :] = np.ma.mean(tasmaxi, axis=0)
                pr_out[i, :, :]     = np.ma.sum(pri, axis=0)
            elif season_type == 'LobellField':
                tas_out[i, :, :]    = LobellFieldGlobalAverage(tasi, area)
                tasmin_out[i, :, :] = LobellFieldGlobalAverage(tasmini, area)
                tasmax_out[i, :, :] = LobellFieldGlobalAverage(tasmaxi, area)
                pr_out[i, :, :]     = LobellFieldGlobalAverage(pri, area)


        prefix = model.lower()
        np.save('growseasons/{0}_{1}_growseason_tas.npy'.format(crop, prefix), tas_out)
        print('tas | saved')

        np.save('growseasons/{0}_{1}_growseason_tasmin.npy'.format(crop, prefix), tasmin_out)
        print('tasmin | saved')

        np.save('growseasons/{0}_{1}_growseason_tasmax.npy'.format(crop, prefix), tasmax_out)
        print('tasmax | saved')

        np.save('growseasons/{0}_{1}_growseason_pr.npy'.format(crop, prefix), pr_out)
        print('pr | saved')
    return True


#######################################################################################################################
#######################################################################################################################
# """
# EXECUTION
# """
#
# tas = loadData('tas')
# tasmin = loadData('tasmin')
# tasmax = loadData('tasmax')
# prec = loadData('pr')
#
# tas     = np.ma.masked_where(tas > 350, tas)
# tas     = np.ma.masked_where(tas < 240, tas)
# tasmin  = np.ma.masked_where(tasmin > 350, tasmin)
# tasmin  = np.ma.masked_where(tasmin < 240, tasmin)
# tasmax  = np.ma.masked_where(tasmax > 350, tasmax)
# tasmax  = np.ma.masked_where(tasmax < 240, tasmax)
# prec    = np.ma.masked_where(prec > 10, prec)
# prec    = np.ma.masked_where(prec < 0, prec)
#
# models = ['LobellField07']
# for model in models:
#     growingSeasonAll(tas, tasmax, tasmin, prec, season_type='LobellField', model=model)
#
# modelstot = ['LPJmL', 'pDSSAT', 'APSIM-UGOE', 'CARAIB', 'EPIC-IIASA', 'EPIC-TAMU', 'GEPIC',
#           'LPJ-GUESS', 'PEPIC', 'PRYSBI2']
#             #LPJ-GUESS only has maize and the two wheats, no soy or rice
# models = ['CARAIB', 'EPIC-IIASA', 'EPIC-TAMU', 'GEPIC', 'PEPIC', 'PRYSBI2']
# for model in models:
#     growingSeasonVar(tas, tasmax, tasmin, prec, season_type='model', model=model)

#######################################################################################################################
#######################################################################################################################

def growingSeasonGDD():
    model = 'pDSSAT'
    prefix = model.lower()
    crop = 'maize'
    netcdf = crop[:3]
    plantdir = '/project2/ggcmi/AgMIP.output/{0}/phase2/{1}/A0/plant-day/'.format(model, crop)
    hardir = '/project2/ggcmi/AgMIP.output/{0}/phase2/{1}/A0/maty-day/'.format(model, crop)
    plant = hfile(
        plantdir + '{0}_agmerra_fullharm_plant-day_{1}_global_annual_1980_2010_C360_T0_W0_N200_A0.nc4'.format(prefix,
                                                                                                              netcdf),
        'r')['plant-day_{0}'.format(netcdf)][:]
    har = \
        hfile(hardir + '{0}_agmerra_fullharm_maty-day_{1}_global_annual_1980_2010_C360_T0_W0_N200_A0.nc4'.format(prefix,
                                                                                                                 netcdf),
              'r')['maty-day_{0}'.format(netcdf)][:]
    print('{0} {1} | loaded'.format(model, crop))

    plant   = np.ma.masked_where(plant == 1e20, plant)
    har     = np.ma.masked_where(har == 1e20, har)
    plant   = np.ma.median(plant, axis=0)
    har     = np.ma.median(har, axis=0)

    pl = (plant).astype(int)
    pl = np.repeat(pl[np.newaxis, :, :], 730, axis=0)
    ha = (har).astype(int)
    ha = np.repeat(ha[np.newaxis, :, :], 730, axis=0)
    ha = pl + ha

    pl = np.where(ha < 365, pl + 365, pl)
    ha = np.where(ha < 365, ha + 365, ha)

    n_xcuts = 4
    x_cuts = np.arange(0, 360, 90)
    dx = 90
    n_ycuts = 4
    y_cuts = np.arange(0, 720, 180)
    dy = 180

    # dayids = np.reshape(np.repeat(np.repeat(np.arange(0, 730, 1), 360, axis=np.newaxis), 720, axis=np.newaxis),
    #                     (730, 360, 720))
    dayids = np.reshape(np.repeat(np.repeat(np.arange(0, 730, 1), dx, axis=np.newaxis), dy, axis=np.newaxis), (730, dx, dy))

    f1 = '/project2/geos39650/ag_data/agmerra/'

    gdd = np.ma.zeros((30, 360, 720))
    startyear = 1981  # this is index 0 in the following loop
    leapyear = 1984  # this is the first leap year, following leaps are +4
    leap = 0  # add a day for each leap year cumulatively through the time loop
    for i in range(1):
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
                tasmini = np.ma.masked_where(tasmini == 1e20, tasmini)
                tasmini = np.ma.harden_mask(tasmini)

                tasmaxi = hfile(f1 + 'tasmax' + '_agmerra_1980-2010.nc4', 'r')['tasmax'][((i) * 365 + leap):((i) * 365 + 730 + leap), xcutj:xendj, ycutk:yendk]
                tasmaxi = np.ma.masked_where(dayids < plj, tasmaxi)
                tasmaxi = np.ma.masked_where(dayids > haj, tasmaxi)
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
                t_hrs = t_hrs - 273.15
                t_hrs[t_hrs > 29] = 29
                t_hrs             = t_hrs - 10
                t_hrs[t_hrs<0]    = 0
                gdd[i, xcutj:xendj, ycutk:yendk] = np.sum(np.sum(t_hrs*(1/24), axis = 3), axis=0)
                gdd[gdd.mask] = 1e20
    np.save('/project2/geos39650/ag_data/gdd/{0}_{1}_gdd.npy'.format(crop, prefix), gdd.data)
    print('gdd | saved')
    return True

growingSeasonGDD()
