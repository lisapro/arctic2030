'''
Created on 31. jan. 2018

@author: ELP

This is a script for making input data for relaxation for the BROM model 
Raw Source data - Geographically sorted data from World Ocean Database 
https://www.nodc.noaa.gov/OC5/WOD/datageo.html
.nc file created in Ocean Data View program

1. Read the nc file, needed variable. 
https://odv.awi.de/
'''

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

def time_to_run():    
    import timeit
    start = timeit.default_timer()
    stop = timeit.default_timer()
    print ('Seconds',stop - start) 
        
def choose_month(ds,m_num,var,clima_var_m,levels,double_int = False,int_num = 1):
    # get 1 month data 
    month_ds = ds.where(ds.date_time.dt.month == m_num, drop=True)
    
    # group by depth, find mean
    try:
        month_mean = month_ds.groupby(month_ds['var1']).mean()     
        month_depth = month_mean.var1.values.copy()
        month_df = month_mean.to_dataframe()
        # check for nans 
        nonnan = len(month_df[var]) - np.count_nonzero(np.isnan(month_df[var])) 
        if nonnan > 10 :       
            # interpolate nans 
            var = (month_df[var].interpolate(kind = 'nearest'))
            if np.count_nonzero(np.isnan(var)) > 0: 
                #find first non nan 
                for n,v in enumerate(var):
                    if np.isnan(var[n]) == False: 
                        pos = n                     
                        for i in np.arange(0,pos):
                            var[i] = var[pos]
                        break                           
            #f = interpolate.UnivariateSpline(month_depth,var,k=3)
            f = interpolate.interp1d(month_depth,var)
            var_m = f(levels)                            
            var_m[var_m < 0] = 0
            
            if double_int == True:
                for n in range(0,int_num):
                    f = interpolate.UnivariateSpline(levels,var_m,k=2)
                    var_m = f(levels)
            print ("month with data num: ", m_num)        
        else: 
            #print ('many nans',m_num,month_df[var])
            var_m = clima_var_m
    except:
        var_m = clima_var_m      
    return var_m,month_ds 

def add_roms_plot(dss,axis,varname):
    funcs = {'Oxygen':'o2', 'Temperature': 'temp',
             'si':'Si','alk': 'Alk',
             'po4':'po4', 'no3':'no3'} 
    
    var_roms = funcs[varname]
    # Filter data: only septembers after 1990
    dss = dss.where( ((dss['time.year']  > 1990) &  (dss['time.month'] == 9)), drop=True)   
    dss['depth'] = dss['depth'].T      
    # Take only september data from ROMS simulation 

    m = dss.groupby(dss.depth).mean()
    axis.scatter(dss[var_roms],dss.depth,  c = '#de7b5c',
            alpha = 0.2, s  = 3, label = 'Sept \nROMS+ERSEM')
    axis.plot(m[var_roms],sorted(dss.depth[0]), 'o', c = '#660033',
            markersize  = 3, label = 'Sept mean \nROMS+ERSEM',zorder = 8 )

def create_relax_array(ncfile,varname,pl,save,levels,axis,int_num = 1,
                        double_int = False,only_clima_mean = False):
    print ('in create relax array')
    funcs = {'Oxygen':'var4', 'Temperature': 'var2',
             'si':'var6','alk': 'var12','chl': 'var10',
             'po4':'var5', 'no3':'var7','pH': 'var9'}  
    var_from_odv = funcs[varname] 
    #print ('varname',varname)      
    # read ncfile as xarray
    ds = xr.open_dataset(ncfile,drop_variables= ('metavar1'))        
    max_depth = np.int(np.max(ds.var1))
    # get only data from 1960 and later
    ds = ds.where(ds.date_time.dt.year > 1950, drop = True) 
    ds = ds.where(ds.var5 < 1.2, drop=True)
    ds = ds.where(ds.var6 < 40, drop=True)        
       
    # group by depth and find mean for each depth 
    clima_mean = ds.groupby(ds['var1']).mean()    
    clima_depth = clima_mean.var1.values.copy()
    clima_df = clima_mean.to_dataframe()
    
    # interpolate nans 
    clima_var = (clima_df[var_from_odv].interpolate(kind = 'nearest'))     
    # interpolations does not work if 1st variable is nan   
    if np.isnan(clima_var[0]) :
        clima_var[0] = clima_var[1]        
    # interpolate values to standard levels    
    f = interpolate.interp1d(clima_depth,clima_var) #,k=3
    #            f = interpolate.interp1d(month_depth,var)
    #        var_m = f(levels)
    clima_var_m = f(levels)
    
    if double_int == True:
        for n in range(0,int_num):
            f1 = interpolate.UnivariateSpline(levels,clima_var_m,k=2)
            clima_var_m = f1(levels)
       
    var_m, var = choose_month(ds,9,var_from_odv,clima_var_m,levels,double_int,int_num)
    var_m[var_m < 0.] = 0 
    #means.append(var_m)
    depths = levels
    clima_means = clima_var_m
    
    if varname == 'Oxygen':    
        clima_var_m = np.array(clima_var_m)*44.6 
        clima_means = np.array(clima_means)*44.6 
        ds[var_from_odv] = ds[var_from_odv]*44.6  
    elif varname == 'alk': 
   
        clima_var_m = np.array(clima_var_m)*1000   
        clima_means = np.array(clima_means)*1000   
        ds[var_from_odv] = ds[var_from_odv]*1000    
          
    if pl == True:
        axis.plot(clima_var_m,levels,'ko--',zorder = 8,
                markersize = 2,label = 'Sept mean WOD')
        axis.scatter(ds[var_from_odv],ds['var1'],alpha = 0.5, 
                c = '#7f7f7f', label = 'Sept WOD', zorder = 7) 

def roms_average_year(var,cmap,vmax,vmin,title) :
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
       
    fig  = plt.figure(figsize=(10,7), dpi=100 )
    fig.suptitle(title)
    gs = gridspec.GridSpec(3,1)
    gs.update(hspace=0.3,top = 0.95,bottom = 0.05, right = 1.02) 
    ax = fig.add_subplot(gs[0]) 
    ax1 = fig.add_subplot(gs[1])
    ax2 = fig.add_subplot(gs[2])   
    ax2.set_title('Mean Temperature year') 
    #ax3 = fig.add_subplot(gs[3])  
       
    ds_1990[var].plot(ax=ax, x = 'time',y = 'depth',cmap=plt.get_cmap(cmap),vmax = vmax,vmin = vmin)
    ds_2002[var].plot(ax=ax1, x = 'time',y = 'depth',cmap=plt.get_cmap(cmap),vmax = vmax,vmin = vmin)  
    cs = ax2.pcolormesh(X,Y,arr,cmap=plt.get_cmap(cmap),vmax = vmax,vmin = vmin) #, cmap=plt.cm.Reds, alpha=1)
    
    ax2.set_ylim(80,0)
    plt.colorbar(cs, ax = ax2)    
    from scipy.interpolate import interp2d

    arr2 = arr.as_matrix(columns = None)
    x = np.array(arr.columns.values) #dayx 
    y = np.array(sorted(d))
    X1,Y1 = np.meshgrid(x,y)    
    
    x_f = np.arange(1,366,1)
    y_f = np.arange(0,81,1)
    X2,Y2 = np.meshgrid(x_f,y_f)
    
    f = interp2d(x,y, arr2, kind ='cubic')
    #cs3 = ax3.pcolormesh(X2,Y2,f(x_f,y_f),cmap=plt.get_cmap(cmap),vmax = 3,vmin = -3) #arr2)
    #ax3.set_ylim(80,0)   
    #plt.colorbar(cs3, ax = ax3)      
      
    ax1.set_ylim(80,0) 
    ax.set_ylim(80,0) 

    
    start = datetime(1989, 1, 1, 23, 59)
    end = datetime(2014, 1, 1, 23, 59)
    rng = pd.date_range(start, end, freq='A-DEC')
    ax.vlines(rng,0,80,linestyle = '--', lw = 0.5)
    ax1.vlines(rng,0,80,linestyle = '--', lw = 0.5)

    plt.show() 
    
def call_arctic() :
    dss = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_each_day.nc')
    levels = sorted(dss.depth.values)
    
    fig  = plt.figure(figsize=(8,8), dpi=100 )
    
    gs = gridspec.GridSpec(2,2)
    gs.update(hspace=0.3,top = 0.95,bottom = 0.05)
    ax = fig.add_subplot(gs[0]) 
    ax1 = fig.add_subplot(gs[1])
    ax2 = fig.add_subplot(gs[2])     
    ax3 = fig.add_subplot(gs[3])     

    ncfile = r'C:/Users/elp/OneDrive/Python_workspace/Relaxation_WOD/src/data_from_WOD_COLLECTION_Laptev.nc'    
    create_relax_array(ncfile,'Oxygen',
                       True, True, levels, ax,1,double_int = True,only_clima_mean = True) # Arctic
    add_roms_plot(dss,ax,'Oxygen')
    ax.set_title(r'O$_2\ \mu M$')  
    
    create_relax_array(ncfile,'po4',
                       True, True, levels,ax1,3, double_int = True,only_clima_mean =True) # Arctic,
    add_roms_plot(dss,ax1,'po4')    
    ax1.set_title(r'PO$_4\ \mu M$') 
        
    create_relax_array(ncfile,'si',
                       True, True, levels,ax2,10, double_int = True,only_clima_mean = True) # Arctic
    add_roms_plot(dss,ax2,'si')    
    ax2.set_title(r'Si $\mu M$') 
 
    create_relax_array(ncfile,'no3',
                       True, True, levels,ax3, 3, double_int = True,only_clima_mean = True) # Arctic
    add_roms_plot(dss,ax3,'no3')       
    ax3.set_title(r'NO$_3\ \mu M$') 

    for axis in (ax,ax1,ax2,ax3):
        axis.set_ylim(90,0)
    l = ax2.legend(facecolor = 'w',framealpha = 0.5)
    l.set_zorder(20)  # put the legend on top
    plt.savefig('Data/WOD_vs_ROMS.png')
    #plt.show( )
 
#call_arctic()   
#roms_average_year('o2')   
roms_average_year('temp','RdBu_r',3,-3,'Temperature C') #'coolwarm')   
