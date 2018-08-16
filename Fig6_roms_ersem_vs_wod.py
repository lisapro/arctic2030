'''
Created on 31. jan. 2018

@author: ELP

This is a script for comparison 
Raw Source data - Geographically sorted data from World Ocean Database 
https://www.nodc.noaa.gov/OC5/WOD/datageo.html
.nc file created in Ocean Data View program

1. Read the nc file, needed variable. 
https://odv.awi.de/
'''

import seaborn as sns
import os 
from statsmodels.discrete.tests.test_constrained import junk
from matplotlib import gridspec as gs
from scipy.interpolate import UnivariateSpline ,griddata 
from scipy import interpolate 
from netCDF4 import Dataset,num2date, date2num
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
import numpy as np       
from datetime import datetime,time  
import xarray as xr
import numpy.ma as ma
import pandas as pd 

sns.set()
        
def choose_month(ds,m_num,var,clima_var_m,levels,double_int = False,int_num = 1):
    # get 1 month data 
    month_ds = ds.where(ds.date_time.dt.month == m_num, drop=True)
    
    try: # group by depth, find mean
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

            f = interpolate.interp1d(month_depth,var)
            var_m = f(levels)                            
            var_m[var_m < 0] = 0
            
            if double_int == True:
                for n in range(0,int_num):
                    f = interpolate.UnivariateSpline(levels,var_m,k=2)
                    var_m = f(levels)      
        else: 
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
    dss = dss.where((dss['time.month'] == 9) |  (dss['time.month'] == 10), 
                    drop=True) 
    dss['depth'] = dss['depth'].T      
    # Take only september data from ROMS simulation 
    m = dss.groupby(dss.depth).mean()
    for n in range(0,len(dss.time),1):
        axis.scatter(dss[var_roms][n],dss.depth,  c = 'y', #de7b5c',
            alpha = 0.2, s  = 3, label = 'Sept \nROMS+ERSEM')
    axis.plot(m[var_roms],sorted(dss.depth), 'o', c = '#660033',
            markersize  = 3, label = 'Sept mean \nROMS+ERSEM',zorder = 8 )
    
def add_brom_plot(dss,axis,varname):
    funcs = {'Oxygen':'o2', 'Temperature': 'temp',
             'si':'Si', #'alk': 'Alk',
             'po4':'PO4', 'no3':'NO3'} 
    dss['depth'] = dss['depth'].T  #  | (dss['time.month'] == 10)
    dss = dss.where(((dss['time.month'] == 9) ), drop=True)      
    var_brom = funcs[varname]   

    for n in range(0,len(dss.time),2):
        axis.scatter(dss[var_brom][n],dss.depth,  c ='#de7b5c', # '#4f542a',
            alpha = 0.5, s  = size-20, zorder = 9) #, label = 'Sept \nBROM')
         
    # Take only september data from BROM simulation
    m = dss.groupby(dss.depth).mean()    
    axis.plot(m[var_brom],sorted(dss.depth), 'o--', c = '#6d3c2c',
            markersize  = 2, label = 'Sept mean \nBROM',zorder = 9)
         

def add_amk(axis,var):
    funcs = {'Oxygen': 'O2ml', 'Temperature': False,
             'si':'Si', #'alk': 'Alk',
             'po4':'PO4', 'no3':'NO3'}     
    file = os.path.join('Data\seep_data_amk.txt') 
    df = pd.read_csv(file,delimiter = '\t')    
    df['O2ml'] *= 44.6 
    v = funcs[var]
    if v == False: 
        pass 
    else: 
        axis.scatter(df[v],df.depth,c = 'k',s = size,zorder = 10)
      
def get_data_wod(ncfile,varname,pl,save,levels,axis,int_num = 1,
                        double_int = False,only_clima_mean = False):
    funcs = {'Oxygen':'var4', 'Temperature': 'var2',
             'si':'var6','alk': 'var12','chl': 'var10',
             'po4':'var5', 'no3':'var7','pH': 'var9'}  
    var_from_odv = funcs[varname]      
    # read ncfile as xarray
    ds = xr.open_dataset(ncfile,drop_variables= ('metavar1'))        
    max_depth = np.int(np.max(ds.var1))
    
    # get only data from 1950 and later
    ds = ds.where(ds.date_time.dt.year > 1940, drop = True) 
    
    # remove all stations (for all variables) 
    # where silicates > 40 micromoles    
    ds = ds.where(ds.var6 < 15, drop=True)   
    #ds = ds.where(ds.var1 < 81, drop=True)    
    ds = ds.where(ds.var5 < 2, drop=True)
    ds = ds.where(ds.var2 < 3, drop=True)  
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
    f = interpolate.interp1d(clima_depth,clima_var) 
    clima_var_m = f(levels)
    
    if double_int == True:
        for n in range(0,int_num):
            f1 = interpolate.UnivariateSpline(levels,clima_var_m,k=2)
            clima_var_m = f1(levels)
       
    var_m, var = choose_month(ds,9,var_from_odv,
                              clima_var_m,levels,
                              double_int,int_num)
    var_m[var_m < 0.] = 0 
    depths = levels
    clima_means = clima_var_m
    
    if varname == 'Oxygen':    
        clima_var_m = np.array(clima_var_m)*44.6 
        #clima_means = np.array(clima_means)*44.6 
        ds[var_from_odv] = ds[var_from_odv]*44.6  
        
    elif varname == 'alk':   
        clima_var_m = np.array(clima_var_m)*1000   
        #clima_means = np.array(clima_means)*1000   
        ds[var_from_odv] = ds[var_from_odv]*1000    
    
    return  clima_var_m,levels,ds[var_from_odv],ds['var1']    


def plot_data_wod(ncfile,varname,pl,save,levels,axis,int_num = 1,
                        double_int = False,only_clima_mean = False,
                        plt_mean=False): 

    clima_var_m,levels,var_from_odv,depth_raw = get_data_wod(ncfile,varname,pl,
                                                save,levels,axis,int_num,
                                                double_int,only_clima_mean)
    if plt_mean:
        axis.plot(clima_var_m,levels,'ko--',zorder = 8,
            markersize = 2,label = 'Sept mean WOD')
    axis.scatter(var_from_odv,depth_raw,alpha = 0.5, s = size,
            c = '#5b5b5b', label = 'Sept WOD', zorder = 9)
     
def plt_brom_ersem_wod(save = False) :
    
    dss = xr.open_dataset('Data\Laptev_baseline.nc')
    dss_brom_init = xr.open_dataset('Data\water.nc')
    
    levels = sorted(dss.depth.values)
    fig  = plt.figure(figsize=(7,7), dpi=100 )
    
    '''
    Data from World Ocean Database 
    https://www.nodc.noaa.gov/OC5/WOD/datageo.html 
    '''
    ncfile = (r'Data\Laptev_WOD_area_70_150_f_deeper.nc')   
    vars = ['Oxygen','po4','si','no3']
    titles = ['O$_2','PO$_4','Si $','NO$_3']
        
    axes = []
    cols = 2
    gs = gridspec.GridSpec(len(vars) // cols , cols)
    gs.update(hspace=0.3,wspace= 0.3,top = 0.92,bottom = 0.02,right = 0.97)
    
    x_text = 0.0
    y_text = 1.1
    labels = ('A) ','B) ','C) ', 'D) ')
    interps = [1,3,10,3] #levels of interpolation     
    for n in range(0,4):
        row = (n // cols)
        col = n % cols    
        axes.append(fig.add_subplot(gs[row, col]))        
        add_brom_plot(dss,axes[n],vars[n])
        add_amk(axes[n],vars[n])
        #add_roms_plot(dss_roms,axes[n],vars[n])     
        plot_data_wod(ncfile,vars[n],
                True, True, levels, axes[n],
                interps[n], double_int = True,
                only_clima_mean = True)
        axes[n].set_ylim(100,0)
        axes[n].xaxis.set_ticks_position('top')
        axes[n].set_ylabel('depth, m')
        axes[n].set_title(r'{}\ \mu M$'.format(titles[n]),y = 1.1) 
        axes[n].text(x_text, y_text, labels[n], transform=axes[n].transAxes,
                 fontsize=14, fontweight='bold', va='top', ha='right')
    if save == True: 
        plt.savefig('Data/Figure6_WOD_vs_BROM.png')
    else:    
        plt.show()
  
size = 30
    
if __name__ == '__main__':    
    plt_brom_ersem_wod(save = True)   
