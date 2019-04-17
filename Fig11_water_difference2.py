'''
Created on 28. jun. 2017

@author: ELP
'''

import os,sys,datetime 
#from PyQt5 import QtWidgets,QtGui, QtCore
#from PyQt5.QtWidgets import QTableWidget,QTableWidgetItem
from netCDF4 import Dataset,num2date,date2num,date2index

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy.ma as ma

import matplotlib.dates as mdates
#import tkinter as tk 
 
import matplotlib.ticker as ticker
import seaborn as sns
import matplotlib.colors as colors
#from matplotlib.mlab import bivariate_normal
from dateutil.relativedelta import relativedelta
#from bokeh.colors import color
sns.set() 
sns.axes_style({ 'axes.grid': True,'axes.spines.right': True})
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

def read_var(name):

    var_water = np.array(
        fh_water.variables[name][:]).T  
    var_water_base  = np.array(
        fh_water_base.variables[name][:]).T     
    var_water_dif = var_water-var_water_base              
    data_units = fh_water.variables[name].units                  
    return var_water_dif,data_units

from matplotlib.ticker import FormatStrFormatter  
def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
         

def make_plot(water_fname,water_base_fname,case,save = True):      
    figure = plt.figure(figsize=(7.27, 8.1), dpi=100,
                            facecolor='None',edgecolor='None')   

    global fh_water,fh_water_base
    fh_water_base = Dataset(water_base_fname)

    fh_water =  Dataset(water_fname)   

    depth_water = np.array(fh_water.variables['z'][:])        
    min_water = np.amin(depth_water)
    max_water = np.amax(depth_water)
            
    #ice_time_format  = num2date(ice_time) 

    time = fh_water.variables['time']      
    time2 = fh_water.variables['time'][:]
    time_units = fh_water.variables['time'].units
    format_time = num2date(time2,units = time_units,calendar= 'standard')

    start_year = 1991
    stop_year = 1992 

    to_start = datetime.datetime(start_year,1,1,12,0)
    to_stop= datetime.datetime(stop_year,1,1,12,0)
            
    start = date2index(to_start,time,
                        calendar=None, select='nearest')
    stop = date2index(to_stop, time,
                        calendar=None, select='nearest') 
    fontsize = 14
    gs0 = gridspec.GridSpec(4, 1)
    gs0.update(left=0.1, right= 0.9,top = 0.96,bottom = 0.03,
                    wspace=0.03, hspace=0.35)  
    w = 0.01
    h = 0.
    gs = gridspec.GridSpecFromSubplotSpec(2, 2, wspace=w,hspace=h,
        subplot_spec=gs0[0],width_ratios=[17, 1], height_ratios=[1, 4])
    gs1 = gridspec.GridSpecFromSubplotSpec(2, 2, wspace=w,hspace=h,
        subplot_spec=gs0[1],width_ratios=[17, 1], height_ratios=[1, 4])
    gs2 = gridspec.GridSpecFromSubplotSpec(2, 2, wspace=w,hspace=h,
        subplot_spec=gs0[2],width_ratios=[17, 1], height_ratios=[1, 4])
    gs3  = gridspec.GridSpecFromSubplotSpec(2, 2, wspace=w,hspace=h,
        subplot_spec=gs0[3],width_ratios=[17, 1], height_ratios=[1, 4])
    
    #add subplots
    ax0 = figure.add_subplot(gs[0])   
    ax1 = figure.add_subplot(gs[2]) 
    cbax1 = figure.add_subplot(gs[:0, -1]) 

    ax0_a = figure.add_subplot(gs1[0])   
    ax1_a = figure.add_subplot(gs1[2]) 
    cbax1_a = figure.add_subplot(gs1[:0, -1])

    ax0_b = figure.add_subplot(gs2[0])   
    ax1_b = figure.add_subplot(gs2[2]) 
    cbax1_b = figure.add_subplot(gs2[:0, -1])

    ax0_c = figure.add_subplot(gs3[0])  
    ax1_c = figure.add_subplot(gs3[2]) 
    cbax1_c = figure.add_subplot(gs3[:0, -1])
                
    ice_data = Dataset('Data\Laptev_average_year_3year.nc')    
    ice = ice_data.variables['hice'][:]        
    ice_time = ice_data.variables['time']   
        
    start_ice = date2index(to_start, ice_time,
                        calendar=None, select='nearest')
    stop_ice = date2index(to_stop, ice_time,
                        calendar=None, select='nearest') 
        
    units = ice_data.variables['time'].units
    form_ice_time = num2date(ice_time[start_ice:stop_ice],units = units)
    X_water,Y_water = np.meshgrid(time2[start:stop],depth_water)
    X_water = num2date(X_water,units = time_units)  

    ice_data.close()    

    for axis in (ax0,ax0_a,ax0_b,ax0_c):
        CS1 = axis.stackplot(form_ice_time,ice[start_ice:stop_ice],color = '#7B98A8')
        axis.set_xlim(to_start,to_stop)   
        axis.set_ylim(0,max(ice[start_ice:stop_ice]))    
        axis.set_yticklabels([])
        axis.set_xticklabels([]) 
        axis.yaxis.set_label_coords(-0.06, 0.5)
        axis.set_ylabel('Ice', fontsize = fontsize) 

    def add_colorbar(CS,axis):# ,ma1
        if axis == cbax1_c:
            cb = plt.colorbar(CS,cax = axis,
            format = FormatStrFormatter('%.3f'))
        elif axis in (cbax1_a,cbax1_b):   
            cb = plt.colorbar(CS,cax = axis,      
            format = FormatStrFormatter('%.1f')) 
        elif axis == cbax1:   
            cb = plt.colorbar(CS,cax = axis,      
            format = FormatStrFormatter('%i'))        
        return cb

    def to_plot(vname,axis,cb_axis,cmap_water):
        var_water,data_units = read_var(vname)

        mm = abs(np.max(var_water[:,start:stop]))
        mm2 = abs(np.min(var_water[:,start:stop]))
        mm_tot = round(max(mm,mm2),2)
        
        #bounds = np.linspace(-mm_tot, mm_tot, 40)     
        #norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)        
            
        #CS4 = axis.pcolor(X_water,Y_water,var_water[:,start:stop],#norm = norm, #vmin=-mm_tot, vmax=mm_tot,
        #                    cmap = cmap_water) 
      
        axis.contour(X_water,Y_water,var_water[:,start:stop],10,linewidths=0.1, colors='k')
        CS4 = axis.contourf(X_water,Y_water,var_water[:,start:stop],100,#norm = norm, #vmin=-mm_tot, vmax=mm_tot,
                            cmap = cmap_water) 


        ma1 = ma.max(var_water[:,start:stop])
        cb1 = add_colorbar(CS4,cb_axis) #,ma1
        axis.set_ylim(max_water,min_water)  
        axis.xaxis.set_major_formatter(
            mdates.DateFormatter('%b')) 
        axis.yaxis.set_label_coords(-0.06, 0.5)
        axis.set_ylabel('Depth,m', fontsize = fontsize) 

    ccmap_water = plt.get_cmap('RdBu_r')  
    import cmocean
    to_plot('B_CH4_CH4',ax1,cbax1,sns.cubehelix_palette(as_cmap=True,light = 1,dark = 0)) #
    to_plot('B_BIO_O2_rel_sat',ax1_a,cbax1_a,plt.get_cmap('Blues_r')) 
    to_plot('B_C_DIC',ax1_b,cbax1_b,plt.get_cmap('pink_r')) 
    to_plot('B_pH_pH',ax1_c,cbax1_c,sns.cubehelix_palette(start=2,as_cmap=True,light = 0,dark = 1)) 

    ax0.set_title(  r'$A)\ scenario\ {}\ , \Delta\ CH_4\ \mu  M $'.format(case),fontsize = fontsize) 
    ax0_a.set_title(r'$B)\ scenario\ {}\ , \Delta\  O_2\ saturation\ \%  $'.format(case),fontsize = fontsize)
    ax0_b.set_title(r'$C)\ scenario\ {}\ , \Delta\  DIC\ \mu  M $'.format(case),fontsize = fontsize) 
    ax0_c.set_title(r'$D)\ scenario\ {}\ , \Delta\   pH$'.format(case),fontsize = fontsize)   

    fh_water.close()
    
    if save == True: 
        plt.savefig('Fig/Figure11_{}.png'.format(case), format='png') #, dpi=300, transparent=True)             
    elif save == False:      
        plt.show()

base = r'C:\Users\elp\OneDrive - NIVA\BROM_linux_output'
#base_nocap = r'{}\no_capping'.format(base)
base_cap = r'{}\with_capping'.format(base)

water_base_b_fname = r'{}\Baseline_B\water.nc'.format(base)
water_base_o_fname = r'{}\Baseline_O\water.nc'.format(base)

water_B  = r'{}\{}\water.nc'.format(base_cap,'B-Basic-seep')      
water_FR = r'{}\{}\water.nc'.format(base_cap,'FR-Reduced-flux')
water_MI = r'{}\{}\water.nc'.format(base_cap,'MI-Increased-horizontal-mixing')       
water_MR = r'{}\{}\water.nc'.format(base_cap,'MR-Reduced-horizontal-mixing')
water_OI = r'{}\{}\water.nc'.format(base_cap,'OI-Increased-Oxidation-Rate')
water_S  = r'{}\{}\water.nc'.format(base_cap,'S-Small-bubbles')      

if __name__ == '__main__':

    make_plot(water_B, water_base_b_fname,case = 'B',save = False)      
    '''make_plot(water_FR,water_base_b_fname,case = 'FR')       
    make_plot(water_MR,water_base_b_fname,case = 'MR')  
    make_plot(water_MI,water_base_b_fname,case = 'MI')     
    make_plot(water_S, water_base_b_fname,case = 'S')          
    make_plot(water_OI,water_base_o_fname,case = 'OI')'''