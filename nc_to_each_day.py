'''
Created on 12. sep. 2017

@author: ELP
'''
# script read netcdf file created in main.py script
# write and interpolate the  data to each day 

import os, time
import datetime
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d,interp2d

import numpy as np

f = Dataset('ROMS_Laptev_Sea_NETCDF3_CLASSIC.nc')

# dimensions
ocean_time =  f.variables['time'][:]
depth =  f.variables['depth'][:]
depth2 =  f.variables['depth2'][:]

# 2d
sal = f.variables['sal'][:][:]
temp = f.variables['temp'][:][:]
kz =  f.variables['Kz_s'][:][:]
rho =  f.variables['rho'][:][:]

# 1d
hice =  f.variables['hice'][:]
snow_thick =  f.variables['snow_thick'][:]
tisrf =  f.variables['tisrf'][:]

times = num2date(ocean_time,units= 'seconds since 1948-01-01 00:00:00' ,
                                   calendar= 'standard')

newtimes =  np.arange(ocean_time[0], ocean_time[-1:], 86400)
newdates = num2date(newtimes,units= 'seconds since 1948-01-01 00:00:00',
                                   calendar= 'standard')

#interpolate 1D variables 
to_interp_hice = interp1d(ocean_time, hice)
to_interp_snow_thick = interp1d(ocean_time, snow_thick)
to_interp_tisrf = interp1d(ocean_time, tisrf)
    
interp_hice = to_interp_hice(newtimes)
interp_snow_thick = to_interp_snow_thick(newtimes)
interp_tisrf = to_interp_tisrf(newtimes)

# plot to test 1d
def plot_1d():
    plt.plot(newdates,interp_hice,'o',alpha = 0.7, label = 'interp')
    #plt.plot(times,hice,'--',alpha = 0.9, label = 'raw')
    plt.legend()
    plt.show()

    



# interpolate 2D variables 
to_interp_sal = interp2d(ocean_time,depth,sal.T,kind = 'linear')
interp_sal = to_interp_sal(newtimes,depth)
interp_sal = interp_sal.T

to_interp_temp = interp2d(ocean_time,depth,temp.T,kind = 'linear')
interp_temp = to_interp_temp(newtimes,depth)
interp_temp = interp_temp.T

to_interp_rho = interp2d(ocean_time,depth,rho.T,kind = 'linear')
interp_rho = to_interp_rho(newtimes,depth)
interp_rho = interp_rho.T

to_interp_kz = interp2d(ocean_time,depth2,kz.T,kind = 'linear')
interp_kz = to_interp_kz(newtimes,depth2)
interp_kz = interp_kz.T


#plot to test 2d
Times,V_depth = np.meshgrid(ocean_time,depth)
Times2,V_depth2 = np.meshgrid(newtimes,depth)
Dates2 = num2date(Times2,units= 'seconds since 1948-01-01 00:00:00',
                                   calendar= 'standard') 
Dates = num2date(Times,units= 'seconds since 1948-01-01 00:00:00',
                                   calendar= 'standard') 
plt.pcolormesh(newtimes,depth2,interp_temp.T) #,
#               vmin=0, vmax=0.1) #, alpha = 0.9, label = 'interp')

#plt.legend()Dates2,V_depth2,c = 
plt.show()




#Write data to new netcdf
#### Create nectdf file 
def write_nc():
    # coordinates of needed station 
    st_lon = 126.82
    st_lat = 76.77
    nc_format = 'NETCDF3_CLASSIC' 
    f1 = Dataset('ROMS_Laptev_Sea_{}_each_day.nc'.format(nc_format), mode='w', format= nc_format)
    
    f1.description= (
        "lat=%3.2f,lon=%3.2f file from ROMS data interpolated to 1day timedelta"%(st_lat,st_lon))
    
    f1.source = 'Elizaveta Protsenko (elp@niva.no)'
    
    f1.history = 'Created ' + time.ctime(time.time())
    
    f1.createDimension('time', len(newtimes))
    f1.createDimension('z', len(depth))
    f1.createDimension('z2', len(depth2))
    
    
    v_depth2 = f1.createVariable('depth2','f8',('z2',), zlib= False)
    v_depth2.long_name = "Z-depth matrix for kz, direction down" ;
    v_depth2.units = "meter"
    v_depth2[:]= depth2
    
    v_depth = f1.createVariable('depth','f8',('z',), zlib= False)
    v_depth.long_name = "Z-depth matrix, direction down" ;
    v_depth.units = "meter"
    v_depth[:]= depth
    
    
    
    v_time = f1.createVariable('time', 'f8', ('time',), zlib=False)
    v_time.long_name = 'Time in seconds since 1948-01-01 00:00:00'
    v_time.units = 'seconds since 1948-01-01 00:00:00'
    v_time.field = 'time, scalar, series'
    v_time.calendar='standard'
    v_time[:] = newtimes
    
    v_temp=f1.createVariable('temp', 'f8', ('time','z'), zlib=False)
    v_temp.long_name = "Ocean temperature"
    v_temp.units = "degree Celsius"
    v_temp[:,:] = interp_temp
    
    v_sal = f1.createVariable('sal', 'f8', ('time','z'), zlib=False)
    v_sal.long_name = "Time-averaged salinity"
    v_sal.units = "psu"
    v_sal[:,:] = interp_sal
    
    v_kz = f1.createVariable('Kz_s', 'float', ('time','z2'), zlib=False)
    v_kz.long_name = 'Salinity vertical diffusion coefficient'
    v_kz.units = 'meter2 second-1'
    #v_kz.coordinates = ("time, depth2")
    v_kz[:,:] = interp_kz #.T
    
    v_rho = f1.createVariable('rho', 'f8', ('time','z'), zlib=False)
    v_rho.long_name = 'time-averaged density anomaly'
    v_rho.units = 'kilogram meter-1'
    v_rho[:,:] = interp_rho
    
    v_hice = f1.createVariable('hice', 'f8', ('time',), zlib=False)
    v_hice.long_name = 'time-averaged ice thickness in cell'
    v_hice.units = 'meter'
    v_hice[:] = interp_hice
    
    v_snow_thick = f1.createVariable('snow_thick', 'f8', ('time',), zlib=False)
    v_snow_thick.long_name = 'time-averaged thickness of snow cover'
    v_snow_thick.units = 'meter'
    v_snow_thick[:] = interp_hice
    
    v_tisrf = f1.createVariable('tisrf', 'f8', ('time',), zlib=False)
    v_tisrf.long_name = 'time-averaged temperature of ice surface'
    v_tisrf.units = 'degree Celsius'
    v_tisrf[:] = interp_tisrf
    
    f1.close()

#write_nc()