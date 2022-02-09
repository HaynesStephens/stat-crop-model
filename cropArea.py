from h5py import File as hfile

def getarea(crop):
    # CAL = calorie per mass (not area)
    # nvar = crop
    # rm = rainfed area
    # rmI = irrigated area
    filedir = '/project2/ggcmi/AgMIP.output/Jim_Emulator/agmerra/'
    if crop   == 'winter_wheat':
        netcdf = 'wwh'
        nvar   = 'wheat'
        rm     = hfile(filedir + 'winter_and_spring_wheat_areas_v1_180627.nc4')['wwh_rf_area'][:]
        rmI    = hfile(filedir + 'winter_and_spring_wheat_areas_v1_180627.nc4')['wwh_ir_area'][:]
        CAL    = 334 * 0.87
    elif crop == 'spring_wheat':
        netcdf = 'swh'
        nvar   = 'wheat'
        rm     = hfile(filedir + 'winter_and_spring_wheat_areas_v1_180627.nc4')['swh_rf_area'][:]
        rmI    = hfile(filedir + 'winter_and_spring_wheat_areas_v1_180627.nc4')['swh_ir_area'][:]
        CAL    = 334 * 0.87
    else:
        netcdf = crop[:3]
        nvar   = crop
        rm     = hfile(filedir + '%s.nc4'%(crop))['rainfed'][:]
        rmI    = hfile(filedir + '%s.nc4'%(crop))['irrigated'][:]
        if crop  == 'rice':   CAL = 280 * 0.88
        elif crop == 'maize': CAL = 356 * 0.88
        else:                 CAL = 335 * 0.91

    rm[rm   > 1000000] = 0
    rmI[rmI > 1000000] = 0
    rm[rm   < 1] = 0
    rmI[rmI < 1] = 0
    return(rm, rmI, CAL, netcdf, nvar)
