from growingSeasonVar import growingSeasonVar
from cropArea import getarea
import numpy as np


def LobellFieldGrowingSeason(crop):
    print('Lobell&Field {0} | loading'.format(crop))
    pl = np.zeros((360, 720))
    ha = np.zeros((360, 720))
    if crop == 'maize':
        pl = pl + 181
        ha = ha + 242
    elif crop == 'soy':
        pl = pl + 181
        ha = ha + 242
    elif crop == 'rice':
        pl = pl + 0
        ha = ha + 303
    elif crop == 'winter_wheat' or crop == 'spring_wheat':
        pl = pl + 120
        ha = ha + 303
    return pl, ha


def LobellFieldGlobalAverage(arr_i, area):
    weights = area.copy()
    arr_mean = np.ma.mean(arr_i, axis=0)
    weights[arr_mean == 0] = 0
    weights[arr_mean.mask] = 0
    return np.ma.sum(arr_mean * weights) / np.ma.sum(weights)


if __name__ == '__main__':
    crops = ['maize', 'soy', 'rice', 'winter_wheat', 'spring_wheat']
    variables = ['tasmin', 'tasmax', 'pr']
    for crop in crops:
        if crop == 'winter_wheat':
            netcdf = 'wwh'
        elif crop == 'spring_wheat':
            netcdf = 'swh'
        else:
            netcdf = crop[:3]

        for var in variables:
            pl, ha = LobellFieldGrowingSeason(crop)
            rm, rmI, CAL, netcdf, nvar = getarea(crop)
            area = rm + rmI
            growingSeasonVar(var, lambda x: LobellFieldGlobalAverage(x, area), pl, ha, 'lobell2007_{0}_{1}'.format(netcdf, var))



