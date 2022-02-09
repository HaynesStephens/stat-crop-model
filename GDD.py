#!/bin/env python
from helper import *
import os
import fnmatch
from h5py import File as hfile
import numpy as np
import pandas as pd
import scipy.odr as odr


# see if mins > maxes
# just swap mix and max if min>max
# get GDD through Roberts method
tmax_hr = np.reshape(np.repeat(tmax, 24, axis=np.newaxis), (730, 360, 720, 24))
tmin_hr = np.reshape(np.repeat(tmin, 24, axis=np.newaxis), (730, 360, 720, 24))
hrs     = np.reshape(np.tile(np.arange(24), (730, 360, 720)), (730, 360, 720, 24))
cos_hrs = np.cos(hrs*np.pi/12)
t_hrs   = ((-1) * ((tmax_hr - tmin_hr)/2) * cos_hrs) + tmin_hr





###########################################################################################################################

np.random.seed(1234)

temp = hfile('agmerra/dailydata/tas_agmerra_1980-2010.nc4', 'r')['tas'][:]


def main(crop, temp):
    models = ['LPJmL']

    for model in models:
        if crop == 'winter_wheat':
            netcdf = 'wwh'
        elif crop == 'spring_wheat':
            netcdf = 'swh'
        else:
            netcdf = crop[:3]
        if model == 'JULES':
            ag = 'wfdei';   nit = 'NNA'
        elif model == 'CARAIB':
            ag = 'agmerra'; nit = 'NNA'
        elif model == 'PROMET':
            ag = 'erai';    nit = 'N200'
        else:
            ag = 'agmerra'; nit = 'N200'

        rootdir = '/project2/ggcmi/AgMIP.output/%s/phase2/%s/A0/yield/' % (model, crop)
        plantdir = '/project2/ggcmi/AgMIP.output/%s/phase2/%s/A0/plant-day/' % (model, crop)
        hardir = '/project2/ggcmi/AgMIP.output/%s/phase2/%s/A0/maty-day/' % (model, crop)

        filelist = os.listdir(rootdir)
        files = fnmatch.filter(filelist, '*C360*W0*%s*.nc4' % (nit))

        ids = np.reshape(np.repeat(np.repeat(np.arange(0, 730, 1), 360, axis=np.newaxis), 720, axis=np.newaxis),
                         (730, 360, 720))
        # f1    = '/project/ggcmi/AgMIP.output/LPJmL/phase2/fasttrack/HadGEM2-ES/rcp8p5/maize/c360_n200/'
        # plant = hfile( f1 + 'lpjml_HadGEM2-ES_rcp8p5_fullharm_plant-day_%s_global_annual_1951_2099_noirr2.nc4'%( netcdf ), 'r' )[ 'plant-day_%s'%( netcdf ) ][:]
        # har   = hfile( f1 + 'lpjml_HadGEM2-ES_rcp8p5_fullharm_maty-day_%s_global_annual_1951_2099_noirr2.nc4'%( netcdf ),  'r' )[ 'maty-day_%s'%( netcdf ) ][:]
        # plant = np.nan_to_num( plant )
        # har   = np.nan_to_num( har )
        # plant[plant > 365] = 0
        # har[har     > 365] = 0
        # plant = plant[31:, :, :]
        # har   = har[31:, :, :]

        plantdir = '/project2/ggcmi/AgMIP.output/LPJmL/phase2/%s/A0/plant-day/' % (crop)
        hardir = '/project2/ggcmi/AgMIP.output/LPJmL/phase2/%s/A0/maty-day/' % (crop)

        plant = \
        hfile(plantdir + 'lpjml_agmerra_fullharm_plant-day_mai_global_annual_1980_2010_C360_T0_W0_N200_A0.nc4', 'r')[
            'plant-day_%s' % (netcdf)][:]
        har = hfile(hardir + 'lpjml_agmerra_fullharm_maty-day_mai_global_annual_1980_2010_C360_T0_W0_N200_A0.nc4', 'r')[
                  'maty-day_%s' % (netcdf)][:]

        plant = np.nan_to_num(plant)
        har = np.nan_to_num(har)
        plant[plant > 365] = 0
        har[har > 365] = 0
        plant = np.ma.median(plant, axis=0)
        har = np.ma.median(har, axis=0)

        pl = (plant).astype(int)
        pl = np.repeat(pl[np.newaxis, :, :], 730, axis=0)
        ha = (har).astype(int)
        ha = np.repeat(ha[np.newaxis, :, :], 730, axis=0)
        ha = pl + ha
        pl = np.where(ha < 365, pl + 365, pl)
        ha = np.where(ha < 365, ha + 365, ha)

        for file in files:
            I = file.split('T')[1]
            I = I.split('_')[0]

            hdd = np.zeros((30, 360, 720))
            gdd = np.zeros((30, 360, 720))
            T = np.zeros((30, 360, 720))

            for i in range(30):
                tas = temp[(i * 365):(i * 365 + 730), :, :] + int(I)
                tas = np.ma.masked_where(ids < pl, tas)
                tas = np.ma.masked_where(ids > ha, tas)
                tas = tas - 273.15
                tas = np.ma.masked_where(tas < -100, tas)
                tasH = np.ma.masked_where(tas < 29, tas)
                tasH = tasH - 29

                tasG = np.ma.masked_where(tas < 10, tas)
                tasG = np.ma.masked_where(tasG > 29, tasG)
                tasG = tasG - 10

                T[i, :, :] = np.ma.mean(tas, axis=0)
                hdd[i, :, :] = np.ma.sum(tasH, axis=0)
                gdd[i, :, :] = np.ma.sum(tasG, axis=0)

            np.save('%s_%s_%s_phaseII_HDD_1982GS' % (crop, model, I), hdd)
            np.save('%s_%s_%s_phaseII_GDD_1982GS' % (crop, model, I), gdd)
            np.save('%s_%s_%s_phaseII_T_phaseII_GS' % (crop, model, I), T)


main('maize', temp)
