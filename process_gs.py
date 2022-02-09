import numpy as np
from growingSeasonVar import growingSeasonVar


def modelGrowingSeason(model, crop):
    prefix = model.lower()
    if crop == 'winter_wheat':
        netcdf = 'wwh'
    elif crop == 'spring_wheat':
        netcdf = 'swh'
    else:
        netcdf = crop[:3]
    plantdir = '/project2/ggcmi/AgMIP.output/{0}/phase2/{1}/A0/plant-day/'.format(model, crop)
    hardir = '/project2/ggcmi/AgMIP.output/{0}/phase2/{1}/A0/maty-day/'.format(model, crop)
    print('{0} {1} | loading'.format(model, crop))
    plant = hfile(
        plantdir + '{0}_agmerra_fullharm_plant-day_{1}_global_annual_1980_2010_C360_T0_W0_N200_A0.nc4'.format(prefix,
                                                                                                              netcdf),
        'r')['plant-day_{0}'.format(netcdf)][:]
    har = \
        hfile(hardir + '{0}_agmerra_fullharm_maty-day_{1}_global_annual_1980_2010_C360_T0_W0_N200_A0.nc4'.format(prefix,
                                                                                                                 netcdf),
              'r')['maty-day_{0}'.format(netcdf)][:]
    print('{0} {1} | loaded'.format(model, crop))

    plant = np.nan_to_num(plant)
    har = np.nan_to_num(har)
    plant[plant > 365] = 0
    har[har > 365] = 0
    plant = np.ma.median(plant, axis=0)
    har = np.ma.median(har, axis=0)
    pl = (plant).astype(int)
    ha = (har).astype(int)
    ha = pl + ha

    print('{0} {1} | seasons determined'.format(model, crop))
    return pl, ha


if __name__ == '__main__':
    models= ['pDSSAT', 'PEPIC', 'GEPIC', 'EPIC-TAMU', 'LPJmL', 'LPJ-GUESS', 'CARAIB', 'PROMET', 'JULES']
    crops = ['maize', 'soy', 'rice', 'winter_wheat', 'spring_wheat']
    variables = ['tas', 'tasmin', 'tasmax', 'pr']
    for model in models:
        for crop in crops:
            if crop == 'winter_wheat':
                netcdf = 'wwh'
            elif crop == 'spring_wheat':
                netcdf = 'swh'
            else:
                netcdf = crop[:3]

            for var in variables:
                pl, ha = modelGrowingSeason(model, crop)
                if var == 'pr':
                    aggfxn = lambda x: np.ma.sum(x, axis=0)
                else:
                    aggfxn = lambda x: np.ma.mean(x, axis=0)
                growingSeasonVar(var, aggfxn, pl, ha, '{0}_{1}_{2}'.format(model, netcdf, var))

