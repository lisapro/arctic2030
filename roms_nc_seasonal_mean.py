import os 
from statsmodels.discrete.tests.test_constrained import junk
from matplotlib import gridspec as gs
from scipy.interpolate import UnivariateSpline  
from scipy import interpolate 
from netCDF4 import Dataset,num2date, date2num
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
import numpy as np       
import datetime
from datetime import time  
#import pandas as pd
import xarray as xr
from scipy import interpolate 
import numpy.ma as ma
from scipy.interpolate import griddata
import pandas as pd 
import os, time
from statsmodels.nonparametric.smoothers_lowess import lowess

    
def get_averaged_value(var):
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True)         
    #d = ds.depth.values    
    arr = ds[var].to_pandas().sort_index() 
    arr = arr.resample('D').ffill()    
    arr['dayofyear'] = arr.index.dayofyear  
    arr = arr.pivot_table(arr,index = 'dayofyear', aggfunc=np.median).T.sort_index() 
    arr = arr.reindex(index=arr.index[::-1])
    arr = arr.loc[:,1:365] 
    return arr.T


def get_averaged_value_2_year(var):
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True)         
    #d = ds.depth.values    
    arr = ds[var].to_pandas().sort_index() 
    arr = arr.resample('D').ffill()    
    arr['dayofyear'] = arr.index.dayofyear  
    arr2 = arr.copy()
    arr2['dayofyear'] = arr2['dayofyear'] + 365
        
    arr = arr.pivot_table(arr,index = 'dayofyear', aggfunc=np.median).T.sort_index() 
    arr2 = arr2.pivot_table(arr2,index = 'dayofyear', aggfunc=np.median).T.sort_index()  
            
    arr = arr.reindex(index=arr.index[::-1])
    arr = arr.loc[:,1:365] 

    arr = pd.concat([arr,arr2],axis = 1)
    arr = arr.loc[:,1:730]    
    
    
   # print (arr.shape)
    return arr.T



def get_averaged_value_1d(var):
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True)        
    ds['dayofyear'] = ds.time.dt.dayofyear
    ds = ds[[var,'dayofyear']]
    arr = ds.to_dataframe()
    #arr = arr.resample('D',how = 'ffill')
    arr = arr.groupby('dayofyear').mean()
    print (arr.shape)
    arr = arr[1:366]

    
    smoothed_arr = lowess(arr[var].values,arr.index.values,frac=0.1)[:,1]
    print(len(smoothed_arr))
    #[:,1] 
    #plt.title('Mean ice thickness,m ')
    #plt.plot(arr.index.values,arr[var].values)
    #plt.plot(arr.index.values,smoothed_arr)
    #plt.show()
    #plt.savefig('meanhice.png', transparent = True)
     

    return smoothed_arr,arr


def get_averaged_value_1d_2_year(var):
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True)        
    ds['dayofyear'] = ds.time.dt.dayofyear
    ds = ds[[var,'dayofyear']]
    arr = ds.to_dataframe()
    #arr = arr.resample('D',how = 'ffill')
    arr = arr.groupby('dayofyear').mean()

    arr = arr[1:366]
    arr2 = arr.append(arr)
    print ('arr2 1d average', len(arr2))
    smoothed_arr = lowess(arr[var].values,arr.index.values,frac=0.05)[:,1]
    if var == 'hice':         
        smoothed_arr_c = []
        for n in smoothed_arr:
            if n < 0.5: 
                n = n - 0.1
            if n < 0:
                n = 0    
            smoothed_arr_c.append(n)   
        smoothed_arr = smoothed_arr_c
        
    smoothed_arr2 = np.append(smoothed_arr,smoothed_arr) 
    
    #[:,1] 
    plt.title(var)
    
    plt.plot(arr.index.values,arr[var].values)
    plt.plot(arr.index.values,smoothed_arr)
    plt.show()
    #plt.savefig('meanhice.png', transparent = True)
     

    return smoothed_arr2,arr2



def get_dimensions():
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True)         
    d = ds.depth.values    
    d2 = ds.depth2.values 
    arr = ds['o2'].to_pandas().sort_index() 
    arr = arr.resample('D').ffill()    
    arr['dayofyear'] = arr.index.dayofyear  
    
    arr = arr.pivot_table(arr,index = 'dayofyear', aggfunc=np.median).T.sort_index() 
    arr = arr.loc[:,1:365] 
    shape = arr.shape
    
    return shape,arr.columns,d,d2

def get_dimensions_2_year():
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True)         
    d = ds.depth.values    
    d2 = ds.depth2.values 
    arr = ds['o2'].to_pandas().sort_index() 
    arr = arr.resample('D').ffill()    
    arr['dayofyear'] = arr.index.dayofyear  
    arr2 = arr.copy()
    arr2['dayofyear'] = arr2['dayofyear'] + 365
    
    arr = arr.pivot_table(arr,index = 'dayofyear', aggfunc=np.median).T.sort_index() 
    arr2 = arr2.pivot_table(arr2,index = 'dayofyear', aggfunc=np.median).T.sort_index()     
   
    arr = arr.loc[:,1:365] 
    arr = pd.concat([arr,arr2],axis = 1)
    #arr = arr.loc[:,1:731]
  
    shape = arr.shape

    return shape,arr.columns,d,d2


def make_nc():
    nc_format = 'NETCDF3_CLASSIC'
    f1 = Dataset('Data\Laptev_average_year.nc', mode='w', format= nc_format)
    f1.description="file from ROMS averaged to one year" 
    f1.source = 'Elizaveta Protsenko (elp@niva.no)'
    f1.history = 'Created ' + time.ctime(time.time())
    
    shape,days, depth,depth2 = get_dimensions()
    
    start = datetime.datetime(1991,1,1)
    format_time = [start]
    seconds = [date2num(start, units = 'seconds since 1948-01-01 00:00:00',calendar = 'standard')]    
    for n in np.arange(0,364):
        delta  = datetime.timedelta(days=int(days[n]))
        date= delta+start
        second = date2num(date, units = 'seconds since 1948-01-01 00:00:00',calendar = 'standard')
        seconds.append(second)
        format_time.append(date)   
    #seconds = date2num(format_time, units = 'seconds since 1948',calendar = 'standard')    
    slb = fig_bub_influx.slblt()
    f1.createDimension('time',  size=len(seconds))
    f1.createDimension('days',  size=shape[1])    
    f1.createDimension('depth', size=shape[0])
    f1.createDimension('depth2', size=len(depth2))
    
    v_depth = f1.createVariable('depth','f8',('depth',), zlib= False)
    v_depth.long_name = "Z-depth matrix, direction up" 
    v_depth.units = "meter"
    v_depth[:] = depth

    v_depth = f1.createVariable('depth2','f8',('depth2',), zlib= False)
    v_depth.long_name = "Z-depth2 matrix, direction up" 
    v_depth.units = "meter"
    v_depth[:] = depth2

      
    v_days = f1.createVariable('days', 'f8', ('days',), zlib=False)
    v_days.long_name = 'number of day in a year'
    v_days.field = 'time, scalar, series'
    v_days[:] = days.values

    v_time = f1.createVariable('time', 'f8', ('days',), zlib=False)
    v_time.long_name = 'Time in seconds since 1948-01-01 00:00:00'
    v_time.units = 'seconds since 1948-01-01 00:00:00'
    v_time.field = 'time, scalar, series'
    v_time.calendar='standard'
    v_time[:] = seconds
       
    v_o2 = f1.createVariable('o2', 'f8', ('time','depth'), zlib=False)
    v_o2.long_name = 'time-averaged oxygen/oxygen '
    v_o2.units = 'mmol O_2/m^3'
    v_o2[:] = get_averaged_value('o2')

    v_temp = f1.createVariable('temp', 'f8', ('time','depth'), zlib=False)
    v_temp.long_name = 'time-averaged ocean temperatue '
    v_temp.units = 'Celsius'
    v_temp[:] = get_averaged_value('temp')

    v_sal = f1.createVariable('sal', 'f8', ('time','depth'), zlib=False)
    v_sal.long_name = 'time-averaged salinity '
    v_sal.units = 'psu'
    v_sal[:] = get_averaged_value('sal')

    v_rho = f1.createVariable('rho', 'f8', ('time','depth'), zlib=False)
    v_rho.long_name = 'time-averaged density anomaly '
    v_rho.units = 'kilogram meter-3'
    v_rho[:] = get_averaged_value('rho')

    v_Kz_s = f1.createVariable('Kz_s', 'f8', ('time','depth2'), zlib=False)
    v_Kz_s.long_name = 'Salinity vertical diffusion coefficient'
    v_Kz_s.units = 'm^2/sec'
    v_Kz_s[:] = get_averaged_value('Kz_s')
    
    v_hice = f1.createVariable('hice', 'f8', ('time',), zlib=False)
    v_hice.long_name = 'time-averaged ice thickness in cell'
    v_hice.units = 'meter'
    v_hice[:] = get_averaged_value_1d('hice')[0]


    v_tisrf = f1.createVariable('tisrf', 'f8', ('time',), zlib=False)
    v_tisrf.long_name = 'time-averaged temperature of ice surface'
    v_tisrf.units = 'Celsius'
    v_tisrf[:] = get_averaged_value_1d('tisrf')[0]

    v_swrad = f1.createVariable('swrad', 'f8', ('time',), zlib=False)
    v_swrad.long_name = 'time-averaged solar shortwave radiation flux'
    v_swrad.units = 'watt meter-2'
    v_swrad.negative_value = 'upward flux,cooling'
    v_swrad[:] = get_averaged_value_1d('swrad')[0]

    v_snow_thick = f1.createVariable('snow_thick', 'f8', ('time',), zlib=False)
    v_snow_thick.long_name = 'time-averaged thickness of snow cover'
    v_snow_thick.units = 'meter'
    v_snow_thick[:] = get_averaged_value_1d('snow_thick')[0]

    v_po4 = f1.createVariable('po4', 'f8', ('time','depth'), zlib=False)
    v_po4.long_name = 'time-averaged phosphate/phosphorus'
    v_po4.units = 'mmol P/m^3'
    v_po4[:] = get_averaged_value('po4')
       
    v_no3 = f1.createVariable('no3', 'f8', ('time','depth'), zlib=False)
    v_no3.long_name = 'time-averaged nitrate/nitrogen'
    v_no3.units = 'mmol N/m^3'
    v_no3[:] = get_averaged_value('no3')
        
    v_Si = f1.createVariable('Si', 'f8', ('time','depth'), zlib=False)
    v_Si.long_name = 'time-averaged silicate/silicate'
    v_Si.units = 'mmol Si/m^3'
    v_Si[:] = get_averaged_value('Si')        
 
    flux_B1,cont_B1 = fig_bub_influx.calculate_scenario_B1(depth,False,days)
    
    v_B1 = f1.createVariable('B1f', 'f8', ('time','depth'), zlib=False)
    v_B1.long_name = 'Methane inflow scenario B1'
    v_B1.units = 'mmol CH4/m^2 sec'
    v_B1[:] = flux_B1 

    v_B1_cont = f1.createVariable('B1c', 'f8', ('time','depth'), zlib=False)
    v_B1_cont.long_name = 'Methane content in bubbles B1'
    v_B1_cont.units = 'mmol CH4 in bubbles'
    v_B1_cont[:] =  cont_B1 

    flux_B1_50,cont_B1_50 = fig_bub_influx.calculate_scenario_B1_50(depth,False,slb,days)
    
    v_B1_50 = f1.createVariable('B1_50f', 'f8', ('time','depth'), zlib=False)
    v_B1_50.long_name = 'Methane inflow scenario B1_50'
    v_B1_50.units = 'mmol CH4/m^2 sec'
    v_B1_50[:] = flux_B1_50 

    v_B1_50_cont = f1.createVariable('B1_50c', 'f8', ('time','depth'), zlib=False)
    v_B1_50_cont.long_name = 'Methane content in bubbles B1'
    v_B1_50_cont.units = 'mmol CH4 in bubbles'
    v_B1_50_cont[:] =  cont_B1_50         

    flux_B1_stop,cont_B1_stop = fig_bub_influx.calculate_scenario_B1_stop(depth,False,slb,days)
    slb_year = fig_bub_influx.calculate_baseline(depth,days)
    
    v_B1_stop = f1.createVariable('B1_stopf', 'f8', ('time','depth'), zlib=False)
    v_B1_stop.long_name = 'Methane inflow scenario B1_stop'
    v_B1_stop.units = 'mmol CH4/m^2 sec'
    v_B1_stop[:] = flux_B1_stop 
 
    v_B1_stop_cont = f1.createVariable('B1_stopc', 'f8', ('time','depth'), zlib=False)
    v_B1_stop_cont.long_name = 'Methane content in bubbles B1_stop'
    v_B1_stop_cont.units = 'mmol CH4 in bubbles'
    v_B1_stop_cont[:] =  cont_B1_stop
 
    v_B1_stop_cont = f1.createVariable('Slb', 'f8', ('time','depth'), zlib=False)
    v_B1_stop_cont.long_name = 'Methane solubility'
    v_B1_stop_cont.units = 'mmol/l CH4 '
    v_B1_stop_cont[:]  = slb_year       
    #ds_1990[var].plot(ax=ax, x = 'time',y = 'depth',cmap=plt.get_cmap(cmap),vmax = vmax,vmin = vmin)
 
    #cs = ax2.pcolormesh(X,Y,arr,cmap=plt.get_cmap(cmap),vmax = vmax,vmin = vmin) 
    
    f1.close()

def make_nc_2_year():
    nc_format = 'NETCDF3_CLASSIC'
    f1 = Dataset('Data\Laptev_average_year_2year.nc', mode='w', format= nc_format)
    f1.description="file from ROMS averaged to one year" 
    f1.source = 'Elizaveta Protsenko (elp@niva.no)'
    f1.history = 'Created ' + time.ctime(time.time())
    
    n_years = 2
    
    shape,days, depth,depth2 = get_dimensions_2_year()

    start = datetime.datetime(1991,1,1)
    format_time = [start]
    seconds = [date2num(start, units = 'seconds since 1948-01-01 00:00:00',calendar = 'standard')]    
    for n in np.arange(0,365*n_years-1):
        delta  = datetime.timedelta(days=int(days[n]))
        date= delta+start
        second = date2num(date, units = 'seconds since 1948-01-01 00:00:00',calendar = 'standard')
        seconds.append(second)
        format_time.append(date)   
        
        
    #print ('seconds',len(seconds))   
    #seconds = date2num(format_time, units = 'seconds since 1948',calendar = 'standard')    
    slb = fig_bub_influx.slblt()
    f1.createDimension('time',  size=len(seconds))
    f1.createDimension('days',  size=shape[1])    
    f1.createDimension('depth', size=shape[0])
    f1.createDimension('depth2', size=len(depth2))
    
    v_depth = f1.createVariable('depth','f8',('depth',), zlib= False)
    v_depth.long_name = "Z-depth matrix, direction up" 
    v_depth.units = "meter"
    v_depth[:] = depth

    v_depth = f1.createVariable('depth2','f8',('depth2',), zlib= False)
    v_depth.long_name = "Z-depth2 matrix, direction up" 
    v_depth.units = "meter"
    v_depth[:] = depth2

      
    v_days = f1.createVariable('days', 'f8', ('days',), zlib=False)
    v_days.long_name = 'number of day in a year'
    v_days.field = 'time, scalar, series'
    v_days[:] = days.values

    v_time = f1.createVariable('time', 'f8', ('time',), zlib=False)
    v_time.long_name = 'Time in seconds since 1948-01-01 00:00:00'
    v_time.units = 'seconds since 1948-01-01 00:00:00'
    v_time.field = 'time, scalar, series'
    v_time.calendar='standard'
    v_time[:] = seconds
       
    v_o2 = f1.createVariable('o2', 'f8', ('time','depth'), zlib=False)
    v_o2.long_name = 'time-averaged oxygen/oxygen '
    v_o2.units = 'mmol O_2/m^3'
    v_o2[:] = get_averaged_value_2_year('o2')

    v_temp = f1.createVariable('temp', 'f8', ('time','depth'), zlib=False)
    v_temp.long_name = 'time-averaged ocean temperatue '
    v_temp.units = 'Celsius'
    v_temp[:] = get_averaged_value_2_year('temp')

    v_sal = f1.createVariable('sal', 'f8', ('time','depth'), zlib=False)
    v_sal.long_name = 'time-averaged salinity '
    v_sal.units = 'psu'
    v_sal[:] = get_averaged_value_2_year('sal')

    v_rho = f1.createVariable('rho', 'f8', ('time','depth'), zlib=False)
    v_rho.long_name = 'time-averaged density anomaly '
    v_rho.units = 'kilogram meter-3'
    v_rho[:] = get_averaged_value_2_year('rho')

    v_Kz_s = f1.createVariable('Kz_s', 'f8', ('time','depth2'), zlib=False)
    v_Kz_s.long_name = 'Salinity vertical diffusion coefficient'
    v_Kz_s.units = 'm^2/sec'
    v_Kz_s[:] = get_averaged_value_2_year('Kz_s')
    
    v_hice = f1.createVariable('hice', 'f8', ('time',), zlib=False)
    v_hice.long_name = 'time-averaged ice thickness in cell'
    v_hice.units = 'meter'
    v_hice[:] = get_averaged_value_1d_2_year('hice')[0]


    v_tisrf = f1.createVariable('tisrf', 'f8', ('time',), zlib=False)
    v_tisrf.long_name = 'time-averaged temperature of ice surface'
    v_tisrf.units = 'Celsius'
    v_tisrf[:] = get_averaged_value_1d_2_year('tisrf')[0]

    v_swrad = f1.createVariable('swrad', 'f8', ('time',), zlib=False)
    v_swrad.long_name = 'time-averaged solar shortwave radiation flux'
    v_swrad.units = 'watt meter-2'
    v_swrad.negative_value = 'upward flux,cooling'
    v_swrad[:] = get_averaged_value_1d_2_year('swrad')[0]

    v_snow_thick = f1.createVariable('snow_thick', 'f8', ('time',), zlib=False)
    v_snow_thick.long_name = 'time-averaged thickness of snow cover'
    v_snow_thick.units = 'meter'
    v_snow_thick[:] = get_averaged_value_1d_2_year('snow_thick')[0]

    v_po4 = f1.createVariable('po4', 'f8', ('time','depth'), zlib=False)
    v_po4.long_name = 'time-averaged phosphate/phosphorus'
    v_po4.units = 'mmol P/m^3'
    v_po4[:] = get_averaged_value_2_year('po4')
       
    v_no3 = f1.createVariable('no3', 'f8', ('time','depth'), zlib=False)
    v_no3.long_name = 'time-averaged nitrate/nitrogen'
    v_no3.units = 'mmol N/m^3'
    v_no3[:] = get_averaged_value_2_year('no3')
        
    v_Si = f1.createVariable('Si', 'f8', ('time','depth'), zlib=False)
    v_Si.long_name = 'time-averaged silicate/silicate'
    v_Si.units = 'mmol Si/m^3'
    v_Si[:] = get_averaged_value_2_year('Si')        
 
    days_1 = np.arange(1,366)
    flux_B1,cont_B1 = fig_bub_influx.calculate_scenario_B1(depth,False,days_1)

    flux_B1 = pd.concat([flux_B1,flux_B1],axis = 0)
    cont_B1 = pd.concat([cont_B1,cont_B1],axis = 0)
    
    v_B1 = f1.createVariable('B1f', 'f8', ('time','depth'), zlib=False)
    v_B1.long_name = 'Methane inflow scenario B1'
    v_B1.units = 'mmol CH4/m^2 sec'
    v_B1[:] = flux_B1 

    v_B1_cont = f1.createVariable('B1c', 'f8', ('time','depth'), zlib=False)
    v_B1_cont.long_name = 'Methane content in bubbles B1'
    v_B1_cont.units = 'mmol CH4 in bubbles'
    v_B1_cont[:] =  cont_B1 

    flux_B1_50,cont_B1_50 = fig_bub_influx.calculate_scenario_B1_50(depth,False,slb,days_1)
    slb_year = fig_bub_influx.calculate_baseline(depth,days_1)
    slb_year = pd.concat([slb_year,slb_year],axis = 0)
    flux_B1_50 = pd.concat([flux_B1_50,flux_B1_50],axis = 0)
    cont_B1_50 = pd.concat([cont_B1_50,cont_B1_50],axis = 0)
        
    v_B1_50 = f1.createVariable('B1_50f', 'f8', ('time','depth'), zlib=False)
    v_B1_50.long_name = 'Methane inflow scenario B1_50'
    v_B1_50.units = 'mmol CH4/m^2 sec'
    v_B1_50[:] = flux_B1_50 

    v_B1_50_cont = f1.createVariable('B1_50c', 'f8', ('time','depth'), zlib=False)
    v_B1_50_cont.long_name = 'Methane content in bubbles B1'
    v_B1_50_cont.units = 'mmol CH4 in bubbles'
    v_B1_50_cont[:] =  cont_B1_50         

    v_B1_slb = f1.createVariable('Slb', 'f8', ('time','depth'), zlib=False)
    v_B1_slb.long_name = 'Methane solubility'
    v_B1_slb.units = 'mmol/l CH4 '
    v_B1_slb[:]  = slb_year      
    
    f1.close()

def plot_roms_average_year(var,title,cmap = 'gist',vmax = None,vmin = None) :
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True)     
    d = ds.depth.values
    
    ds_1990 = ds.where((ds['time.year'] < 2002), drop=True) 
    ds_2002= ds.where((ds['time.year'] >= 2002), drop=True)
    arr = ds[var].to_pandas().sort_index() 
    arr = arr.resample('D').ffill()  
    t = arr.index.dayofyear   
    arr['dayofyear'] = t    
    
    arr = arr.pivot_table(arr,index = 'dayofyear', aggfunc=np.median).T.sort_index() #columns = 'depth', 
    X,Y = np.meshgrid(arr.columns,arr.index)
  
    fig  = plt.figure(figsize=(8,9), dpi=100 )
    fig.suptitle(title)
    gs = gridspec.GridSpec(3,1)
    gs.update(hspace=0.35,top = 0.95,bottom = 0.05, right = 1.05,left = 0.1) 
    ax = fig.add_subplot(gs[0]) 
    ax1 = fig.add_subplot(gs[1])
    ax2 = fig.add_subplot(gs[2])   
    ax2.set_title('Average year') 
       
    ds_1990[var].plot(ax=ax, x = 'time',y = 'depth',cmap=plt.get_cmap(cmap),vmax = vmax,vmin = vmin)
    ds_2002[var].plot(ax=ax1, x = 'time',y = 'depth',cmap=plt.get_cmap(cmap),vmax = vmax,vmin = vmin)  
    cs = ax2.pcolormesh(X,Y,arr,cmap=plt.get_cmap(cmap),vmax = vmax,vmin = vmin) 
    
    plt.colorbar(cs, ax = ax2)    

    for axis in (ax,ax1,ax2):
        axis.set_ylabel('depth, m')
        axis.set_ylim(80,0)
    start = datetime(1989, 1, 1, 23, 59)
    end = datetime(2014, 1, 1, 23, 59)
    rng = pd.date_range(start, end, freq='A-DEC')
    ax.vlines(rng,0,80,linestyle = '--', lw = 0.5)
    ax1.vlines(rng,0,80,linestyle = '--', lw = 0.5)
    ax2.set_xlabel('number of day in a year')
    #plt.savefig(r'C:\Users\elp\OneDrive\Python_workspace\arctic2030\Data\mean_roms_{}.png'.format(var))
    plt.show() 
    plt.clf()


    
'''plot_roms_average_year('temp','Temperature C','RdBu_r',3,-3) #'coolwarm')   
plot_roms_average_year('o2','O$_2\ \mu M$','viridis',400,300)
plot_roms_average_year('po4','PO$_4\ \mu M$','viridis',0.8,0)
plot_roms_average_year('Si','Si $ \mu M$','viridis',15,2)'''
#plot_roms_average_year('no3','NO$_3\ \mu M$','viridis',8,0)
import fig_bub_influx 


#get_averaged_value('o2')
#make_nc()
make_nc_2_year()
