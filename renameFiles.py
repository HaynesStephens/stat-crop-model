import glob
import os



dir = '/project2/geos39650/ag_data/growseasons'
os.chdir(dir)
file_list = glob.glob('*.npy')

# for i in file_list:
#     file_name = i.split('_')
#     if (file_name[0] == 'spring') or (file_name[0] == 'winter'):
#         crop = file_name[0] + '_' + file_name[1]
#         model = file_name[2]
#         gs = file_name[3]
#         var = file_name[4].split('.')[0]
#     else:
#         crop = file_name[0]
#         model = file_name[1]
#         gs = file_name[2]
#         var = file_name[3].split('.')[0]
#     model = model.replace('-', '')
#
#     if crop == 'winter_wheat':
#         netcdf = 'wwh'
#     elif crop == 'spring_wheat':
#         netcdf = 'swh'
#     else:
#         netcdf = crop[:3]
#
#     new_name = '{0}_{1}_{2}.npy'.format(model, netcdf, var)
#
#     os.system('mv {0} {1}'.format(i, new_name))

for i in file_list:
    file_name = i.split('_')
    print(file_name)
    if 'test' not in file_name[0]:
        model = file_name[0]
        netcdf = file_name[1]
        var = file_name[2].split('.')[0]
        if model == 'apsimugoe':
            new_name = '{0}_{1}_{2}.npy'.format('apsim-ugoe', netcdf, var)
            os.system('mv {0} {1}'.format(i, new_name))


