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
from matplotlib import gridspec
sns.set() 

base = r'C:\Users\elp\OneDrive - NIVA\BROM_linux_output'
base_nocap = r'{}\no_capping'.format(base)
base_cap = r'{}\with_capping'.format(base)

def get_path(name):
    return r'{}\{}\water.nc'.format(base_cap,name)  

water_B_cap = get_path('B-Basic-seep')  
water_FR_cap = get_path('FR-Reduced-flux')  
water_MI_cap = get_path('MI-Increased-horizontal-mixing')        
water_MR_cap = get_path('MR-Reduced-horizontal-mixing')       
water_OI_cap = get_path('OI-Increased-Oxidation-Rate') 
water_S_cap = get_path('S-Small-bubbles') 

scenarios = [water_B_cap,water_FR_cap,water_MI_cap,water_MR_cap,water_OI_cap]

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
    
def make_figure_all_scen():
    df_b = xr.open_dataset(scenarios[0])
    pal = sns.cubehelix_palette(as_cmap= True,light = 1,dark= 0)
    z = np.reshape(np.tile(df_b.z.values, 366),(366,45))
    slb = np.array(list(map(lambda t,s,d: calc_methane_saturation(t,s,d), df_b.temp, df_b.sal, z))).T
    X,Y = np.meshgrid(df_b.time.values[:366],df_b.z.values)

    def getvar(path):
       return xr.open_dataset(path)['B_CH4_CH4'].values.T[:,:366] 

    fig = plt.figure(figsize=(6,7)) #subplots(6, 2, sharex=True)
    gs = gridspec.GridSpec(6, 2,width_ratios= (9,1),hspace = 0.4,
                                wspace = 0.1,right = 0.95,
                                left = 0.1,bottom = 0.05,top = 0.97)
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0]) 
    ax3 = plt.subplot(gs[2,0]) 
    ax4 = plt.subplot(gs[3,0]) 
    ax5 = plt.subplot(gs[4,0]) 
    ax6 = plt.subplot(gs[5,0]) 
    ax7 = plt.subplot(gs[:,1]) 

    def perc(var):
        return 100*(var/slb)
    ticks = np.arange(0,100,10) 

    dfs = [getvar(water_B_cap),getvar(water_OI_cap),getvar(water_MI_cap),
           getvar(water_MR_cap),getvar(water_FR_cap),getvar(water_S_cap)]

    CS1 = ax1.contourf(X,Y,perc(dfs[0]),10, vmin = 0, vmax = 100,cmap = pal)
    fig.colorbar(CS1,cax=ax7,ticks=ticks)

    titles = ('A) Scenario B ','B) Scenario OI','C) Scenario MI',
              'D) Scenario MR','E) Scenario FR','G) Scenario S')
    axes = (ax1,ax2,ax3,ax4,ax5,ax6)
    [a.contourf(X,Y,perc(dfs[k]),10, vmin = 0, vmax = 100,cmap = pal) for k,a in enumerate(axes)]
    [a.set_ylim(80,0) for a in axes]
    [a.set_xticks([]) for a in axes[:-1]]
    [a.set_ylabel('Depth,m') for a in axes]
    [a.set_title('{} - ratio to solubility %'.format(titles[n])) for n,a in enumerate(axes)]
    ax6.xaxis.set_major_formatter(
            mdates.DateFormatter('%b')) 
    plt.savefig('Fig/Fig10_to_solubility_perc.png',format = 'png')
    #plt.show()




if __name__ == '__main__':
    #make_figure(water_MI_cap,'MI')
    #make_figure(water_MR_cap,'MR')
    #make_figure(water_B_cap,'B')
    #make_figure(water_F_cap,'F')
    #make_figure(water_OI_cap,'OI')
    #make_figure(water_S_cap,'S')
   make_figure_all_scen()


