

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np       
import xarray as xr
import pandas as pd 
import os, time
import S5_Calculate_scenarios as scen 
from datetime import date,datetime
from netCDF4 import Dataset,num2date, date2num

def cut_roms():   
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_each_day_test.nc')
    ds = ds.where((
        (ds['time.year']  >= 2000)), drop=True)  
    ds = ds.where((
        (ds['time.year']  <= 2010)), drop=True)    
    ds.to_netcdf('Data\ROMS_Laptev_Sea_2000-2010_cut.nc',format = "NETCDF3_CLASSIC" ) 
    
def cut_roms_1y(year,file_path):
   
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_each_day_test.nc')
    ds = ds.where(((ds['time.year']  == year)), drop=True)  
    ds.to_netcdf(file_path,format = "NETCDF3_CLASSIC" )    
    
def repeat_value_N_year(var,n_years):                      
    return  xr.concat([var,var,var],dim='time')   


def make_nc_3_year(year):
    nc_format = 'NETCDF3_CLASSIC'
    f1 = Dataset(
    'Data\Laptev_{}_year_3year.nc'.format(str(yr)), mode='w',
    format= nc_format)
    f1.description=(
    "file from ROMS {} year cut + methane"
    "seeping scenarios from Single Bubble model".format(str(year)))
    f1.source = 'Elizaveta Protsenko (elp@niva.no)'
    f1.history = 'Created ' + time.ctime(time.time())
    
    n_years = 3
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_{}.nc'.format(str(year)))


    # interpolate ice
    from statsmodels.nonparametric.smoothers_lowess import lowess
    import statsmodels
    
    a=lowess(ds['hice'], ds['time'],frac=0.3)
    ds['hice2'] = a[:,1]
    
    depth = ds.depth.values    
    depth2 = ds.depth2.values 
    days = pd.date_range(start = '01/1/2000', periods=365*n_years)
    
    epoch_1948 = datetime(year=1948,month=1,day=1)   
    seconds  = (days - epoch_1948).total_seconds()

    f1.createDimension('time',  size = len(days))
    f1.createDimension('days',  size = len(days))    
    f1.createDimension('depth', size = len(depth))
    f1.createDimension('depth2', size = len(depth2))   

    v_depth = f1.createVariable('depth','f8',('depth',),
                                 zlib= False)
    v_depth.long_name = "Z-depth matrix, direction up" 
    v_depth.units = "meter"
    v_depth[:] = depth

    v_depth = f1.createVariable('depth2','f8',('depth2',),
                                 zlib= False)
    v_depth.long_name = "Z-depth2 (interfaces) matrix, direction up" 
    v_depth.units = "meter"
    v_depth[:] = depth2
      
    v_time = f1.createVariable('time', 'f8', ('time',), 
                               zlib=False)
    v_time.long_name = 'Time in seconds since 1948-01-01 00:00:00'
    v_time.units = 'seconds since 1948-01-01 00:00:00'
    v_time.field = 'time, scalar, series'
    v_time.calendar='standard'
    v_time[:] = seconds.values
    
   
    v_o2 = f1.createVariable('o2', 'f8', ('time','depth'), zlib=False)
    v_o2.long_name = 'time-averaged oxygen/oxygen '
    v_o2.units = 'mmol O_2/m^3'
    v_o2[:] = repeat_value_N_year(ds['o2'],n_years)
    
    v_temp = f1.createVariable('temp', 'f8', ('time','depth'), zlib=False)
    v_temp.long_name = 'time-averaged ocean temperatue '
    v_temp.units = 'Celsius'
    v_temp[:] = repeat_value_N_year(ds['temp'],n_years)

    v_sal = f1.createVariable('sal', 'f8', ('time','depth'), zlib=False)
    v_sal.long_name = 'time-averaged salinity '
    v_sal.units = 'psu'
    v_sal[:] = repeat_value_N_year(ds['sal'],n_years)
    
    v_rho = f1.createVariable('rho', 'f8', ('time','depth'), zlib=False)
    v_rho.long_name = 'time-averaged density anomaly '
    v_rho.units = 'kilogram meter-3'
    v_rho[:] = repeat_value_N_year(ds['rho'],n_years)

    v_Kz_s = f1.createVariable('Kz_s', 'f8', ('time','depth2'), zlib=False)
    v_Kz_s.long_name = 'Salinity vertical diffusion coefficient'
    v_Kz_s.units = 'm^2/sec'
    v_Kz_s[:] = repeat_value_N_year(ds['Kz_s'],n_years)

    v_po4 = f1.createVariable('po4', 'f8', ('time','depth'), zlib=False)
    v_po4.long_name = 'time-averaged phosphate/phosphorus'
    v_po4.units = 'mmol P/m^3'
    v_po4[:] = repeat_value_N_year(ds['po4'],n_years)
       
    v_no3 = f1.createVariable('no3', 'f8', ('time','depth'), zlib=False)
    v_no3.long_name = 'time-averaged nitrate/nitrogen'
    v_no3.units = 'mmol N/m^3'
    v_no3[:] = repeat_value_N_year(ds['no3'],n_years)
        
    v_Si = f1.createVariable('Si', 'f8', ('time','depth'), zlib=False)
    v_Si.long_name = 'time-averaged silicate/silicate'
    v_Si.units = 'mmol Si/m^3'
    v_Si[:] = repeat_value_N_year(ds['Si'],n_years)  
    
    v_hice = f1.createVariable('hice', 'f8', ('time',), zlib=False)
    v_hice.long_name = 'time-averaged ice thickness in cell'
    v_hice.units = 'meter'
    v_hice[:] = repeat_value_N_year(ds['hice2'],n_years)

    
    v_tisrf = f1.createVariable('tisrf', 'f8', ('time',), zlib=False)
    v_tisrf.long_name = 'time-averaged temperature of ice surface'
    v_tisrf.units = 'Celsius'
    v_tisrf[:] =  repeat_value_N_year(ds['tisrf'],n_years)

    v_swrad = f1.createVariable('swrad', 'f8', ('time',), zlib=False)
    v_swrad.long_name = 'time-averaged solar shortwave radiation flux'
    v_swrad.units = 'watt meter-2'
    v_swrad.negative_value = 'upward flux,cooling'
    v_swrad[:] =  repeat_value_N_year(ds['swrad'],n_years)

    v_snow_thick = f1.createVariable('snow_thick', 'f8', ('time',), zlib=False)
    v_snow_thick.long_name = 'time-averaged thickness of snow cover'
    v_snow_thick.units = 'meter'
    v_snow_thick[:] =  repeat_value_N_year(ds['snow_thick'],n_years)

    slb = scen.calculate_baseline(days)/10
    v_B1_slb = f1.createVariable('Slb', 'f8', ('time','depth'), zlib=False)
    v_B1_slb.long_name = 'Methane solubility'
    v_B1_slb.units = 'mmol/l CH4 '
    v_B1_slb[:]  = slb 
                  
    f1.close()
    
def oneyear_3year() :
    yr = 2010
    new_file = 'Data\ROMS_Laptev_Sea_{}.nc'.format(str(yr))
    cut_roms_1y(yr,new_file)     
    make_nc_3_year(yr)


def interpolate_average(): 
    ds = xr.open_dataset('Data\Laptev_average_year_3year.nc')
    # interpolate ice
    from statsmodels.nonparametric.smoothers_lowess import lowess
    import statsmodels
    
    a=lowess(ds['hice'], ds['time'],frac=0.3)
    ds.assign('hice2':a[:,1])    
    #print(ds['hice2'],ds['hice'])
    ds.to_netcdf('Data\Laptev_average_year_3year_av_hice.nc',format = "NETCDF3_CLASSIC" ) 
    
    
interpolate_average()    
    