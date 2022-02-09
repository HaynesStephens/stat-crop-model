### Load plant and harvest days ###
rec = xr.open_dataset(rootdir + 'growing_season_datesv1.25/{0}_rf_growing_season_dates_v1.25.nc4'.format(netcdf))
pla = np.nan_to_num(np.flip(rec['planting day'].values, axis=0), nan=-999)
har = np.nan_to_num(np.flip(rec['harvest day'].values, axis=0), nan=-999)
pla = np.ma.masked_outside(pla, 365, 0).harden_mask()
har = np.ma.masked_outside(har, 365, 0).harden_mask()
