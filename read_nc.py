'''
Created on 2. okt. 2017

@author: ELP
'''
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import matplotlib.pyplot as plt

f = Dataset('ROMS_Laptev_Sea_NETCDF3_CLASSIC_each_day.nc')
depth =  f.variables['depth'][:]
sal = f.variables['sal'][:][7000]
temp = f.variables['temp'][:][7000]
ocean_time =  f.variables['time'][7000]
times = num2date(ocean_time,
                 units= 'seconds since 1948-01-01 00:00:00' ,
                                   calendar= 'standard')

plt.plot(temp,depth)
plt.ylim(150,0)
plt.title(str(times)+', temperature')
plt.show()