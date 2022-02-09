#!/bin/env python
import os; import fnmatch; from netCDF4 import Dataset; import numpy as np; import pandas as pd
def main(model, crop):
    ### write netcdf4 ###
    nc_out       = Dataset( '%s_%s_ggcmi_phase2_emulator_A1.nc4'%(model, crop), 'w', format = 'NETCDF4' )
    nc_out.title = "Aggregated Variable for Regression X"
    nc_out.institution = "University of Chicago"
    nc_out.contact     = "haynes13@uchicago.edu"
    poly_dim = nc_out.createDimension( 'poly', None )
    lat_dim = nc_out.createDimension( 'lat', 360 )
    lon_dim = nc_out.createDimension( 'lon', 720 )
    lat     = nc_out.createVariable( 'lat', 'f8', ('lat',) )
    lat.units     = 'degrees_north'
    lat.long_name = 'latitude'
    lat[:] = np.linspace( 89.75, -89.75, 360 )
    lon    = nc_out.createVariable( 'lon', 'f8', ('lon',) )
    lon.units     = 'NA'
    lon.long_name = 'longitude'
    lon[:] = np.linspace( -179.75, 179.75, 720 )

    y_out  = nc_out.createVariable( 'K_rf', 'f8', ( 'poly', 'lat', 'lon' ), fill_value = -999 )
    y_out.units = 'none'
    y_out.long_name = 'Polynomial parameters for emulation of ggcmi phase II RAINFED yield output (A1): [0, :, :] = intercept, [1, :, :] = K1 ... and so on. See Franke et al. 2020 in Geoscientific Model Development for details.'
    y_out[:, :]     = np.load('A1/%s_%s.npy'%(model, crop))
    y2_out  = nc_out.createVariable( 'K_ir', 'f8', ( 'poly', 'lat', 'lon' ), fill_value = -999 )
    y2_out.units = 'none'
    y2_out.long_name = 'Polynomial parameters for emulation of ggcmi phase II IRRIGATED yield output (A1): [0, :, :] = intercept, [1, :, :] = K1 ... and so on. See Franke et al. 2020 in Geoscientific Model Development for details.'
    y2_out[:, :]     = np.load('A1/%s_%s_I.npy'%(model, crop))
    nc_out.close()
