'''
Created on 28. jun. 2017

@author: ELP
'''
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import os,sys,datetime 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import xarray as xr
import pandas as pd

sns.set() 

base = r'C:\Users\elp\OneDrive - NIVA\BROM_linux_output'
base_nocap = r'{}\no_capping'.format(base)
base_cap = r'{}\with_capping'.format(base)

def calc_methane_saturation(temp,sal,d):
    Kh_theta = 1.4e-5 * 101.325 # mol/m3 Pa #Convert to M/atm  
    coef = 1900 #[K]
    T = temp + 273.15
    k_sal = 10**(-(sal/58)*0.127)
    press = 1 + d/10
    conc_sat = Kh_theta*1.e6*press*np.exp((coef)*(1/T - 1/298.15))*k_sal #Micromol/l
    return conc_sat

def make_figure(scenario,scen_name):
    df_scen = xr.open_dataset(scenario)
   
    z = np.reshape(np.tile(df_scen.z.values, 366),(366,45))
    #print (df_scen.depth2.shape)
    slb = np.array(list(map(lambda t,s,d: calc_methane_saturation(t,s,d), df_scen.temp, df_scen.sal, z))).T
    b = df_scen['B_CH4_CH4'].values.T
    b = b[:,:366] 
    fig, [ax1, ax2, ax3] = plt.subplots(3, 1, sharex=True)
    
    X,Y = np.meshgrid(df_scen.time.values[:366],df_scen.z.values)
    CS = ax1.contourf(X,Y,slb,30)
    fig.colorbar(CS,ax=ax1)
    ax1.set_title('Equlibrium solubility')

    CS2 = ax2.contourf(X,Y,b,30)
    fig.colorbar(CS2,ax=ax2)
    ax2.set_title(r'$CH_4$ concentration scenario {}'.format(scen_name))
     
    perc =  100*(b/slb)
    pal = sns.cubehelix_palette(as_cmap= True,light = 0,dark= 1)
    CS3 = ax3.contourf(X,Y,perc,10, vmin = 0, vmax = 100,cmap = pal)
    fig.colorbar(CS3,ax=ax3,ticks=[0,20,40,60,80,100])
    ax3.set_title('"%" of equilibrium concentration')

    [a.set_ylim(80,0) for a in (ax1,ax2,ax3)]

    plt.savefig('Fig/Fig13_{}'.format(scen_name))
    #plt.show()
    

#water_B_nocap = r'{}\B-Basic-seep\water.nc'.format(base_nocap)  

water_B_cap = r'{}\{}\water.nc'.format(base_cap,'B-Basic-seep')  
water_F_cap = r'{}\{}\water.nc'.format(base_cap,'F-Reduced-flux')  
water_MI_cap = r'{}\{}\water.nc'.format(base_cap,
                'MI-Increased-horizontal-mixing')        
water_MR_cap = r'{}\{}\water.nc'.format(base_cap,
                'MR-Reduced-horizontal-mixing')       

water_OI_cap = r'{}\{}\water.nc'.format(base_cap,
                'OI-Increased-Oxidation-Rate') 

water_S_cap = r'{}\{}\water.nc'.format(base_cap,
                'S-Small-bubbles') 

if __name__ == '__main__':
    #make_figure(water_MI_cap,'MI')
    #make_figure(water_MR_cap,'MR')
    #make_figure(water_B_cap,'B')
    #make_figure(water_F_cap,'F')
    make_figure(water_OI_cap,'OI')
    make_figure(water_S_cap,'S')
   


