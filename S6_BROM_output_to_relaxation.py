from datetime import time  
import os, time
import xarray as xr
from netCDF4 import Dataset,num2date, date2num
import S3_roms_nc_seasonal_mean as make_nc
import S5_Calculate_scenarios as scen 
import numpy as np
import pandas as pd
'''
ds = xr.open_dataset('Data\water_baseline.nc')
depth = ds['middle_layer_depths'].values
depth2 = ds.depth2.values 
Kz_s = ds.Kz_s.values

print (ds)
r'''


''' BROM reads ROMS and adds BBL to it,
in order to read  '''

import datetime
 
start = datetime.datetime(1991,1,1)
sec_start = date2num(start, 
                units = 'seconds since 1948-01-01 00:00:00',
                calendar = 'standard')
def d2n(delta): 
    d  = datetime.timedelta(days=int(delta))  
    s = date2num((d+start), 
                units = 'seconds since 1948-01-01 00:00:00',
                calendar = 'standard')
    return s 





def make_nc_baseline():
    f_brom = Dataset('Data\water_baseline.nc', mode='r')
    f_roms = Dataset('Data\Laptev_average_year_2year.nc', mode='r')
    
    z = f_brom.variables['z'][5:]
    z2 = f_brom.variables['z_faces'][5:]
    
    depth = f_roms.dimensions['depth']
    depth2 = f_roms.dimensions['depth2']
    times = f_roms.dimensions['time']
    seconds = [d2n(delta) for delta in np.arange(1,366*2)]    
    days = [day for day in np.arange(1,366*2)]

    nc_format = 'NETCDF3_CLASSIC'
    f1 = Dataset('Data\Laptev_baseline.nc', mode='w', format= nc_format)
    f1.description="Baseline after BROM run" 
    f1.source = 'Elizaveta Protsenko (elp@niva.no)'
    f1.history = 'Created ' + time.ctime(time.time())
    
    
    f1.createDimension('time',  size=len(times))
    f1.createDimension('days',  size=len(days))    
    f1.createDimension('depth', size= len(depth))
    f1.createDimension('depth2', size=len(depth2))
    
    
    print (len(times),len(seconds),len(days))
    
    v_depth = f1.createVariable('depth','f8',('depth',), zlib= False)
    v_depth.long_name = "Z-depth matrix, direction up" 
    v_depth.units = "meter"
    v_depth[:] = z
    
    v_depth2 = f1.createVariable('depth2','f8',('depth2',), zlib= False)
    v_depth2.long_name = "Z-depth2 matrix, direction up" 
    v_depth2.units = "meter"
    v_depth2[:] = z2

    v_time = f1.createVariable('time', 'f8', ('time',), zlib=False)
    v_time.long_name = 'Time in seconds since 1948-01-01 00:00:00'
    v_time.units = 'seconds since 1948-01-01 00:00:00'
    v_time.field = 'time, scalar, series'
    v_time.calendar='standard'
    v_time[:] = seconds
    
    v_days = f1.createVariable('days', 'f8', ('days',), zlib=False)
    v_days.long_name = 'number of day in a year'
    v_days.field = 'time, scalar, series'
    v_days[:] = days
      
    
    
    
    f1.close()
    f_brom.close()
    f_roms.close()
    #time = f1.dimensions['time']
    #z = f0.variables['depth'][:]
    
    """shape, days, depth,depth2 = make_nc.get_dimensions_2_year()
    days_1 = np.arange(1,367)
    
    if ('swrad' in list(f1.variables)) == False:
        v_swrad = f1.createVariable('swrad', 'f8', ('time',), zlib=False)
        v_swrad.long_name = 'time-averaged solar shortwave radiation flux'
        v_swrad.units = 'watt meter-2'
        v_swrad.negative_value = 'upward flux,cooling'
        v_swrad[:] = make_nc.get_averaged_value_1d_2_year('swrad')[0]
        
    if ('snow_thick' in list(f1.variables)) == False:
        v_snow_thick = f1.createVariable('snow_thick', 'f8', ('time',), zlib=False)
        v_snow_thick.long_name = 'time-averaged thickness of snow cover'
        v_snow_thick.units = 'meter'
        v_snow_thick[:] = make_nc.get_averaged_value_1d_2_year('snow_thick')[0]
    
    if ('B2_10f' in list(f1.variables)) == False:
        flux_B2_10,cont_B1_10 = scen.calculate_scenarios(depth,False,days_1,'B2_10')
        
        flux_B2_10 = pd.concat([flux_B2_10,flux_B2_10],axis = 0).iloc[:731,:] 
        
        v_B2_10 = f1.createVariable('B2_10f', 'f8', ('time','z'), zlib=False)    
        v_B2_10.long_name = 'Methane inflow scenario B2_10'
        v_B2_10.units = 'mmol CH4/m^2 sec'
        v_B2_10[:] = flux_B2_10 """
        
make_nc_baseline()