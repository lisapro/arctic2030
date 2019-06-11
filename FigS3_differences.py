'''
Created on 28. jun. 2017

@author: ELP
'''
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import os,sys,datetime 
from netCDF4 import Dataset,num2date,date2num,date2index
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy.ma as ma
import matplotlib.dates as mdates 
import matplotlib.ticker as ticker
import seaborn as sns
import matplotlib.colors as colors
from dateutil.relativedelta import relativedelta
import xarray as xr
#from bokeh.colors import color
import matplotlib.ticker as mtick
from matplotlib.ticker import MaxNLocator


sns.set() 
sns.set_palette("bright")
sns.set_style("darkgrid")
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
         

base = r'C:\Users\elp\OneDrive - NIVA\BROM_linux_output'
#base_nocap = r'{}\no_capping'.format(base)
base_cap = r'{}\with_capping'.format(base)

def get_path(name):
    return r'{}\{}\water.nc'.format(base_cap,name) 

water_base_b_fname = r'{}\Baseline_B\water.nc'.format(base)
water_base_o_fname = r'{}\Baseline_O\water.nc'.format(base)

water_B  = get_path('B-Basic-seep')      
water_FR = get_path('FR-Reduced-flux')
water_FR2 = get_path('FR2-Reduced-flux')
water_MI = get_path('MI-Increased-horizontal-mixing')       
water_MR = get_path('MR-Reduced-horizontal-mixing')
water_OI = get_path('OI-Increased-Oxidation-Rate')
water_S  = get_path('S-Small-bubbles') 

def plot_scen(name,var,ax1,label,col):
    ds = xr.open_dataset(name)
    ds_base = xr.open_dataset(water_base_b_fname)

    ds = ds.where(ds.time.dt.year < 1992,drop = True)
    ds_base = ds_base.where(ds_base.time.dt.year < 1992,drop = True) 

    ds['dif_var'] = (('time','z'),ds[var].values-ds_base[var].values)
    
    # find max over time   
    group2 = ds['dif_var'].groupby('time')   
    v = group2.max().rolling(time=5, center=True).mean(dim = xr.ALL_DIMS)
    vmin = group2.min().rolling(time=5, center=True).mean(dim = xr.ALL_DIMS) 
    ln = 1 
    ax1.fill_between(group2.groups.keys(),vmin,v,alpha = 0.3,color = col)  
    ax1.plot(group2.groups.keys(),vmin,color = col,linewidth = ln,label = label,)      
    ax1.plot(group2.groups.keys(),v,color = col,linewidth = ln)       
    #ax2.plot(group2.groups.keys(),v,linewidth = 1)


if __name__ == '__main__':


    def call_all_2(var):
        plt.clf()
        fig = plt.figure(figsize=(7.5,4))
        gs = gridspec.GridSpec(2, 1,hspace = 0.2,
                                    wspace = 0.1,right = 0.95,
                                    left = 0.1,bottom = 0.1,top = 0.95)
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1]) 
        #ax3 = plt.subplot(gs[2]) 

        dfs  = {'OI': water_OI,'B':water_B,'FR2':water_FR2,'MR': water_MR,'MI': water_MI,'S': water_S}

        plot_scen(dfs['B'],var,ax1,'B','#57689b')   
        plot_scen(dfs['FR2'],var,ax1,'FR2','#cc584d')
        plot_scen(dfs['OI'],var,ax1,'OI','#008080')  

        plot_scen(dfs['S'],var,ax2,'S','#3f3f3f') 
        plot_scen(dfs['MR'],var,ax2,'MR','#00b159')
        plot_scen(dfs['MI'],var,ax2,'MI','#e5c100')

        ax1.legend(loc = 'best')
        ax2.legend()
        #ax2.xaxis.set_major_locator = mdates.MonthLocator()
        ax2.xaxis.set_major_formatter(
            mdates.DateFormatter('%b'))     
        #ax1.xaxis.set_major_locator = mdates.MonthLocator()
        ax1.xaxis.set_major_formatter(
            mdates.DateFormatter('%b'))         
        ax1.set_title('Daily ranges of differences for {}'.format(var))    

        #plt.show()
        plt.savefig('Fig/Fig_S3_Diff_all_scen_2axis_{}.png'.format(var))

    def call_all_in_one(var):
        
        fig = plt.figure(figsize=(7.8,4))
        gs = gridspec.GridSpec(1, 1,hspace = 0.2,
                                    wspace = 0.1,right = 0.95,
                                    left = 0.1,bottom = 0.1,top = 0.97)
        ax1 = plt.subplot(gs[0])

        colors = ('#008080','#22F148', '#57689b','#cc584d','#e5c100','#3f3f3f')
        names = ('MR','OI','B','FR2','MI','S')
        dfs  = {'OI': water_OI,'B':water_B,'FR2':water_FR2,'MR': water_MR,'MI': water_MI,'S': water_S}

        [plot_scen(dfs[n],var,ax1,n,colors[m]) for m,n in enumerate(names,start = 0)]

        '''plot_scen(water_OI,var,ax1,'OI')   
        plot_scen(water_B,var,ax1,'B')
        plot_scen(water_FR2,var,ax1,'FR2')
        plot_scen(water_MR,var,ax1,'MR')
        plot_scen(water_MI,var,ax1,'MI') 
        plot_scen(water_S,var,ax1,'S') '''   

        ax1.legend(loc = 'best')
        #plt.show()
        plt.savefig('Fig/Fig_S3_Diff_all_scen_{}.png'.format(var))

    #call_all_in_one('B_CH4_CH4')
    vars = ('B_CH4_CH4','B_BIO_O2_rel_sat') #,'B_pH_pH','B_C_CO3','B_BIO_Phy','B_BIO_Het')
    #[call_all_2(var) for var in vars]
    
    call_all_2('B_CH4_CH4')
    call_all_2('B_BIO_O2_rel_sat')
    call_all_2('B_pH_pH')
    call_all_2('B_C_CO3')    
    call_all_2('B_BIO_Phy')      
    call_all_2('B_BIO_Het') 



    '''call_all_in_one('B_BIO_O2_rel_sat')
    call_all_in_one('B_pH_pH')
    call_all_in_one('B_C_CO3')    
    call_all_in_one('B_BIO_Phy')      
    call_all_in_one('B_BIO_Het')  '''

    #print (group)