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
def plt_brom_ersem_wod(save = False) :
    
    dss = xr.open_dataset('Data\Laptev_baseline.nc')
    dss_roms = xr.open_dataset('Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_each_day.nc')
    #print (dss)
    
    levels = sorted(dss.depth.values)
    
    fig  = plt.figure(figsize=(8,8), dpi=100 )
    
    gs = gridspec.GridSpec(2,2)
    gs.update(hspace=0.3,top = 0.95,bottom = 0.05)
    ax = fig.add_subplot(gs[0]) 
    ax1 = fig.add_subplot(gs[1])
    ax2 = fig.add_subplot(gs[2])     
    ax3 = fig.add_subplot(gs[3])   
      
    '''
    Data from World Ocean Database 
    https://www.nodc.noaa.gov/OC5/WOD/datageo.html 
    
    '''
    ncfile = r'C:/Users/elp/OneDrive/Python_workspace/Relaxation_WOD/src/data_from_WOD_COLLECTION_Laptev.nc'   
    vars = ['Oxygen','po4','si','no3']
    axes = [ax,ax1,ax2,ax3]
    for n in range(0,4):
        add_brom_plot(dss,axes[n],var[n])
        add_roms_plot(dss,axes[n],var[n]) 
        
    plot_data_wod(ncfile,'Oxygen',
                       True, True, levels, ax,1,
                       double_int = True,only_clima_mean = True)
    
    #add_roms_plot(dss_roms,ax,'Oxygen')
    #add_brom_plot(dss,ax,'Oxygen')
    ax.set_title(r'O$_2\ \mu M$')  

        
    plot_data_wod(ncfile,'po4',
                       True, True, levels,ax1,3,
                       double_int = True,only_clima_mean =True)
    #add_roms_plot(dss_roms,ax1,'po4')   
    #add_brom_plot(dss,ax1,'po4') 
    ax1.set_title(r'PO$_4\ \mu M$') 
        
    plot_data_wod(ncfile,'si',
                       True, True, levels,ax2,10,
                       double_int = True,only_clima_mean = True)
    add_roms_plot(dss_roms,ax2,'si')    
    add_brom_plot(dss,ax2,'si')
    ax2.set_title(r'Si $\mu M$') 
 
    plot_data_wod(ncfile,'no3',
                       True, True, levels,ax3,3,
                       double_int = True,only_clima_mean = True)
    add_roms_plot(dss_roms,ax3,'no3')  
    add_brom_plot(dss,ax3,'no3')     
    ax3.set_title(r'NO$_3\ \mu M$') 
    
    for axis in (ax,ax1,ax2,ax3):
        axis.set_ylim(90,0)
    l = ax2.legend(facecolor = 'w',framealpha = 0.5)
    l.set_zorder(20)  #put the legend on top

        
    if save == True: 
        plt.savefig('Data/WOD_vs_ROMS.png')
    else:    
        plt.show()
  
def plt_brom_limits():
    dss = xr.open_dataset('Data\water.nc')
    limlight = dss['B_BIO_LimLight']
    for n in range(215,308):
        if (limlight[n][44]) <0.1 :
            print (n)
            plt.plot(limlight[n],dss.z,  c = 'b', #, s  = 3
            alpha = 0.2) #, label = 'Sept \nBROM')    
        #plt.scatter(dss['B_BIO_LimT'][n],dss.z,  c = 'r',
        #    alpha = 0.3, s  = 3) #, label = 'Sept \nBROM')         
        #plt.scatter(dss['B_BIO_LimP'][n],dss.z,  c = 'k',
        #    alpha = 0.3, s  = 3) #, label = 'Sept \nBROM')             
        #plt.scatter(dss['B_BIO_LimN'][n],dss.z,  c = 'k',
        #    alpha = 0.3, s  = 3) #, label = 'Sept \nBROM')      
        #plt.scatter(dss['B_BIO_LimSi'][n],dss.z,  c = 'k',
        #    alpha = 0.3, s  = 3) #, label = 'Sept \nBROM')            
    plt.ylim(80,0)        
    plt.show()
    
if __name__ == '__main__':    
    #plt_brom_ersem_wod()   
    plt_brom_limits()