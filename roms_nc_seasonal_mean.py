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
from datetime import datetime,time  
#import pandas as pd
import xarray as xr
from scipy import interpolate 
import numpy.ma as ma
from scipy.interpolate import griddata
import pandas as pd 
import os, time



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
    
def get_averaged_value(var):
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True)         
    #d = ds.depth.values    
    arr = ds[var].to_pandas().sort_index() 
    arr = arr.resample('D').ffill()    
    arr['dayofyear'] = arr.index.dayofyear  
    arr = arr.pivot_table(arr,index = 'dayofyear', aggfunc=np.median).T.sort_index() 
    arr = arr.reindex(index=arr.index[::-1])
    ###arr = arr.loc[:,1:365] 
    return arr.T

def get_dimensions():
    ds = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc')
    ds = ds.where((ds['time.year']  >= 1989), drop=True)         
    d = ds.depth.values    
    arr = ds['o2'].to_pandas().sort_index() 
    arr = arr.resample('D').ffill()    
    arr['dayofyear'] = arr.index.dayofyear  
    
    arr = arr.pivot_table(arr,index = 'dayofyear', aggfunc=np.median).T.sort_index() 
    #####arr = arr.loc[:,1:365] 
    shape = arr.shape
    return shape,arr.columns,d

def make_nc():
    nc_format = 'NETCDF3_CLASSIC'
    f1 = Dataset('Data\Laptev_average_year.nc', mode='w', format= nc_format)
    f1.description="file from ROMS averaged to one year" 
    f1.source = 'Elizaveta Protsenko (elp@niva.no)'
    f1.history = 'Created ' + time.ctime(time.time())
    
    shape,days, depth = get_dimensions()
    slb = fig_bub_influx.slblt()
    f1.createDimension('time',  size=shape[1])
    f1.createDimension('depth', size=shape[0])
    
    v_depth = f1.createVariable('depth','f8',('depth',), zlib= False)
    v_depth.long_name = "Z-depth matrix, direction up" 
    v_depth.units = "meter"
    v_depth[:] = depth
      
    v_time = f1.createVariable('time', 'f8', ('time',), zlib=False)
    v_time.long_name = 'number of day in a year'
    v_time.field = 'time, scalar, series'
    v_time[:] = days.values
       
    v_o2 = f1.createVariable('o2', 'f8', ('time','depth'), zlib=False)
    v_o2.long_name = 'time-averaged oxygen/oxygen '
    v_o2.units = 'mmol O_2/m^3'
    v_o2[:] = get_averaged_value('o2')

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
 
    flux_B1,cont_B1 = fig_bub_influx.calculate_scenario_B1(depth,False)
    
    v_B1 = f1.createVariable('B1f', 'f8', ('time','depth'), zlib=False)
    v_B1.long_name = 'Methane inflow scenario B1'
    v_B1.units = 'mmol CH4/m^2 sec'
    v_B1[:] = flux_B1 

    v_B1_cont = f1.createVariable('B1c', 'f8', ('time','depth'), zlib=False)
    v_B1_cont.long_name = 'Methane content in bubbles B1'
    v_B1_cont.units = 'mmol CH4 in bubbles'
    v_B1_cont[:] =  cont_B1 

    flux_B1_50,cont_B1_50 = fig_bub_influx.calculate_scenario_B1_50(depth,False,slb)
    
    v_B1_50 = f1.createVariable('B1_50f', 'f8', ('time','depth'), zlib=False)
    v_B1_50.long_name = 'Methane inflow scenario B1_50'
    v_B1_50.units = 'mmol CH4/m^2 sec'
    v_B1_50[:] = flux_B1_50 

    v_B1_50_cont = f1.createVariable('B1_50c', 'f8', ('time','depth'), zlib=False)
    v_B1_50_cont.long_name = 'Methane content in bubbles B1'
    v_B1_50_cont.units = 'mmol CH4 in bubbles'
    v_B1_50_cont[:] =  cont_B1_50         


    flux_B1_stop,cont_B1_stop = fig_bub_influx.calculate_scenario_B1_stop(depth,False,slb)
    slb_year = fig_bub_influx.calculate_baseline(depth)
    
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
    
'''plot_roms_average_year('temp','Temperature C','RdBu_r',3,-3) #'coolwarm')   
plot_roms_average_year('o2','O$_2\ \mu M$','viridis',400,300)
plot_roms_average_year('po4','PO$_4\ \mu M$','viridis',0.8,0)
plot_roms_average_year('Si','Si $ \mu M$','viridis',15,2)'''
#plot_roms_average_year('no3','NO$_3\ \mu M$','viridis',8,0)
import fig_bub_influx 

#get_averaged_value('o2')
make_nc()