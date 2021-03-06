import os 
#from statsmodels.discrete.tests.test_constrained import junk
from matplotlib import gridspec as gs
from scipy.interpolate import UnivariateSpline  
from scipy import interpolate 
from netCDF4 import Dataset,num2date, date2num
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#from mpl_toolkits.basemap import Basemap
import numpy as np       
import datetime
#import pandas as pd
import xarray as xr
from scipy import interpolate 
import numpy.ma as ma
from scipy.interpolate import griddata
import pandas as pd 
import os, time
from statsmodels.nonparametric.smoothers_lowess import lowess
#import x_fig_bub_influx as fig_bub_influx 
import S5_Calculate_scenarios as scen 
    
def get_averaged_value(var):
    ''' average roms 20 years to one year '''
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
    ''' average roms 20 years to one year,repeat 2 times '''
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
            
    
    arr = arr.loc[:,1:365] 
    arr = pd.concat([arr,arr2],axis = 1)
    arr = arr.loc[:,1:731]    
    arr = arr.reindex(index=arr.index[::-1])
    
   # print (arr.shape)
    return arr.T

def get_averaged_value_N_year(var,n_years):
    ''' average roms over 20 years to one year,repeat N times '''
    if n_years == 1: 
        var =  get_averaged_value(var)    
    elif n_years == 2: 
        var =  get_averaged_value_2_year(var)
    else:             
        ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
        ds = ds.where((ds['time.year']  >= 1989), drop=True)          
        arr = ds[var].to_pandas().sort_index() 
        arr = arr.resample('D').ffill()  
        idx = pd.date_range('1989-01-01 ', '2014-12-31')   
        arr = arr.reindex(idx,method='ffill')
        arr = arr.interpolate()
        arr['dayofyear'] = arr.index.dayofyear    
    
        arr = arr.pivot_table(arr,index = 'dayofyear', aggfunc=np.median).T.sort_index()  
        arr = arr.loc[:,1:365]        
        arr1 = arr.copy() 
        
        for n in range(1,n_years):
            arr1.columns = arr1.columns+365            
            arr = pd.concat([arr,arr1],axis = 1,join='outer')      

        arr = arr.reindex(index=arr.index[::-1]).T     
    return  arr     

def get_averaged_value_1d(var):
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True)        
    ds['dayofyear'] = ds.time.dt.dayofyear
    ds = ds[[var,'dayofyear']]
    arr = ds.to_dataframe()
    arr = arr.groupby('dayofyear').mean()
    arr = arr[1:366]
    smoothed_arr = lowess(arr[var].values,arr.index.values,frac=0.1)[:,1]
    return smoothed_arr,arr

def get_averaged_value_1d_2_year(var):
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True)        
    ds['dayofyear'] = ds.time.dt.dayofyear
    ds = ds[[var,'dayofyear']]
    arr = ds.to_dataframe().groupby('dayofyear').mean()
    arr = arr
    arr2 = arr.append(arr)
    arr2 = arr2[:731]
    smoothed_arr = lowess(arr[var].values,arr.index.values,frac=0.05)[:,1]
    if var == 'hice':         
        smoothed_arr_c = []
        for n in smoothed_arr:
            #if n < 0.5: 
            n = n - 0.15
            if n < 0:
                n = 0    
            smoothed_arr_c.append(n)   
        smoothed_arr = smoothed_arr_c
        
    smoothed_arr2 = np.append(smoothed_arr,smoothed_arr) 
    smoothed_arr2 = smoothed_arr2[:731]
    return smoothed_arr2,arr2

def get_averaged_value_1d_3_year(var):
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True)        
    ds['dayofyear'] = ds.time.dt.dayofyear
    ds = ds[[var,'dayofyear']]
    arr = ds.to_dataframe().groupby('dayofyear').mean()
    smoothed_arr = lowess(arr[var].values,arr.index.values,frac=0.05)[:,1]
    if var == 'hice':         
        smoothed_arr_c = []
        for n in smoothed_arr:
            #if n < 0.5: 
            n = n - 0.15
            if n < 0:
                n = 0    
            smoothed_arr_c.append(n)   
        smoothed_arr = smoothed_arr_c
        
    smoothed_arr2 = np.append(smoothed_arr,smoothed_arr) 
    smoothed_arr2 = np.append(smoothed_arr2,smoothed_arr) 
    smoothed_arr2 = smoothed_arr2[:365*3]
    return smoothed_arr2 

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
    shape = arr.shape

    return shape,arr.columns,d,d2


def get_dimensions_N_year(n_years):
            
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True) 
    d = ds.depth.values    
    d2 = ds.depth2.values          
    arr = ds['o2'].to_pandas().sort_index() 
    arr = arr.resample('D').ffill()  
    idx = pd.date_range('1989-01-01 ', '2014-12-31')   
    arr = arr.reindex(idx,method='ffill').interpolate()
    arr['dayofyear'] = arr.index.dayofyear    

    arr = arr.pivot_table(arr,index = 'dayofyear', aggfunc=np.median).T.sort_index()  
    arr = arr.loc[:,1:365]    
    arr1 = arr.copy() 
    
    for n in range(1,n_years):
        arr1.columns = arr1.columns+365            
        arr = pd.concat([arr,arr1],axis = 1,join='outer')             
    arr = arr.reindex(index=arr.index[::-1])            
    return  arr.shape,arr.columns,d,d2

def make_nc_2_year():
    nc_format = 'NETCDF3_CLASSIC'
    f1 = Dataset(
    'Data\Laptev_average_year_3year.nc', mode='w', format= nc_format)
    f1.description=(
    "file from ROMS averaged + methane seeping scenarios from Single Bubble model")
    f1.source = 'Elizaveta Protsenko (elp@niva.no)'
    f1.history = 'Created ' + time.ctime(time.time())
    
    n_years = 3
    shape,days, depth,depth2 = get_dimensions_N_year(n_years)

    start = datetime.datetime(1991,1,1)
    format_time = [start]
    
    seconds = [date2num(start, 
            units = 'seconds since 1948-01-01 00:00:00',calendar = 'standard')]  
      
    for n in np.arange(1,365*n_years):
        delta  = datetime.timedelta(days=int(days[n]))
        date= delta+start
        second = date2num(date,
                units = 'seconds since 1948-01-01 00:00:00',calendar = 'standard')
        seconds.append(second)
        format_time.append(date)   
        
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
    v_o2[:] = get_averaged_value_N_year('o2',n_years)

    v_temp = f1.createVariable('temp', 'f8', ('time','depth'), zlib=False)
    v_temp.long_name = 'time-averaged ocean temperatue '
    v_temp.units = 'Celsius'
    v_temp[:] = get_averaged_value_N_year('temp',n_years)

    v_sal = f1.createVariable('sal', 'f8', ('time','depth'), zlib=False)
    v_sal.long_name = 'time-averaged salinity '
    v_sal.units = 'psu'
    v_sal[:] = get_averaged_value_N_year('sal',n_years)

    v_rho = f1.createVariable('rho', 'f8', ('time','depth'), zlib=False)
    v_rho.long_name = 'time-averaged density anomaly '
    v_rho.units = 'kilogram meter-3'
    v_rho[:] = get_averaged_value_N_year('rho',n_years)

    v_Kz_s = f1.createVariable('Kz_s', 'f8', ('time','depth2'), zlib=False)
    v_Kz_s.long_name = 'Salinity vertical diffusion coefficient'
    v_Kz_s.units = 'm^2/sec'
    v_Kz_s[:] = get_averaged_value_N_year('Kz_s',n_years)
    
    v_hice = f1.createVariable('hice', 'f8', ('time',), zlib=False)
    v_hice.long_name = 'time-averaged ice thickness in cell'
    v_hice.units = 'meter'
    v_hice[:] = get_averaged_value_1d_3_year('hice')


    v_tisrf = f1.createVariable('tisrf', 'f8', ('time',), zlib=False)
    v_tisrf.long_name = 'time-averaged temperature of ice surface'
    v_tisrf.units = 'Celsius'
    v_tisrf[:] = get_averaged_value_1d_3_year('tisrf')

    v_swrad = f1.createVariable('swrad', 'f8', ('time',), zlib=False)
    v_swrad.long_name = 'time-averaged solar shortwave radiation flux'
    v_swrad.units = 'watt meter-2'
    v_swrad.negative_value = 'upward flux,cooling'
    v_swrad[:] = get_averaged_value_1d_3_year('swrad')

    v_snow_thick = f1.createVariable('snow_thick', 'f8', ('time',), zlib=False)
    v_snow_thick.long_name = 'time-averaged thickness of snow cover'
    v_snow_thick.units = 'meter'
    v_snow_thick[:] = get_averaged_value_1d_3_year('snow_thick')

    v_po4 = f1.createVariable('po4', 'f8', ('time','depth'), zlib=False)
    v_po4.long_name = 'time-averaged phosphate/phosphorus'
    v_po4.units = 'mmol P/m^3'
    v_po4[:] = get_averaged_value_N_year('po4',n_years)
       
    v_no3 = f1.createVariable('no3', 'f8', ('time','depth'), zlib=False)
    v_no3.long_name = 'time-averaged nitrate/nitrogen'
    v_no3.units = 'mmol N/m^3'
    v_no3[:] = get_averaged_value_N_year('no3',n_years)
        
    v_Si = f1.createVariable('Si', 'f8', ('time','depth'), zlib=False)
    v_Si.long_name = 'time-averaged silicate/silicate'
    v_Si.units = 'mmol Si/m^3'
    v_Si[:] = get_averaged_value_N_year('Si',n_years)        
 
    days_1 = np.arange(1,366)
    
    stop = 365*n_years
    zeros = scen.calculate_spin_up(depth,days_1)
       
    '''flux_B1,cont_B1 = scen.calculate_scenarios(depth,days_1,'B1') 
    flux_B1 = (pd.concat([zeros,flux_B1,flux_B1],
            axis = 0,ignore_index=True)).iloc[:stop,:]
    cont_B1 = pd.concat([zeros,cont_B1,cont_B1],axis = 0).iloc[:stop,:]'''

    '''slb_year = scen.calculate_baseline(depth)/10
    slb_year = pd.concat([slb_year,slb_year,slb_year],axis = 0).iloc[:stop,:]
    v_B1_slb = f1.createVariable('Slb', 'f8', ('time','depth'), zlib=False)
    v_B1_slb.long_name = 'Methane solubility to atmosphere '
    v_B1_slb.units = 'mmol/l CH4 '
    v_B1_slb[:]  = slb_year  '''
        
    f1.close()

def plot_roms_average_year(var,title,cmap = 'gist',vmax = None,vmin = None) :
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True)     
    d = ds.depth.values
    
    ds_1990 = ds.where((ds['time.year'] < 1992), drop=True) 
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
    start = datetime.datetime(1989, 1, 1, 23, 59)
    end = datetime.datetime(2014, 1, 1, 23, 59)
    rng = pd.date_range(start, end, freq='A-DEC')
    ax.vlines(rng,0,80,linestyle = '--', lw = 0.5)
    ax1.vlines(rng,0,80,linestyle = '--', lw = 0.5)
    ax2.set_xlabel('number of day in a year')
    #plt.savefig(r'C:\Users\elp\OneDrive\Python_workspace\arctic2030\Data\mean_roms_{}.png'.format(var))
    plt.show() 
    plt.clf()


if __name__ == '__main__':
        
    #get_averaged_value('o2')
    #make_nc()
    make_nc_2_year()
    #plot_roms_average_year('temp','Temperature C','RdBu_r',3,-3) #'coolwarm')   
    '''plot_roms_average_year('po4','PO$_4\ \mu M$','viridis',0.8,0)
    plot_roms_average_year('Si','Si $ \mu M$','viridis',15,2)'''
    #plot_roms_average_year('no3','NO$_3\ \mu M$','viridis',8,0)
    ##plot_roms_average_year('o2','O$_2\ \mu M$','viridis',400,300)
    
    