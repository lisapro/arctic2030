import pandas as pd
import matplotlib.pyplot as plt
from glob import glob 
import numpy as np 
import matplotlib.gridspec as gridspec
from scipy import interpolate 
import xarray as xr
from statsmodels.nonparametric.smoothers_lowess import lowess
import statsmodels
import S1_methane_equilibrium  as me 
import seaborn as sns
sns.set()
global bub_path 
bub_path = (
r"C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re(1)\all_for_sbm_79_mm_tab.dat")

def plot_rad_evol(save = False):
    fig = plt.figure(figsize =(6,4), dpi=150)
    gs = gridspec.GridSpec(2,1)
    gs.update(top = 0.95,bottom = 0.15,hspace = 0.4,right = 0.99)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    df = pd.read_csv(bub_path,delimiter = '\t' )       
    gb = df.groupby(df.radius)
    step = 0.25    
    radii = np.arange(0.5,2,step)
    radii2 = np.arange(2,6.75,step)
    radii3 = np.arange(0.5,6.75,step)
    
    perc_to_atm,perc_to_wat  = [],[]
    maxs,mins = [],[]
    
    for rad in radii: 
        group =  gb.get_group(rad)[['met_cont','depth','rad_evol']] 
        f  = interpolate.interp1d(group.depth,group.rad_evol)
        new_depth = np.arange(group.depth.min(),group.depth.max(),3)
        rs = f(new_depth) 
        s_rad = (rs/np.max(rs))*100 
        
        ax0.plot(group.rad_evol,group.depth,alpha = 0.2,
                 linewidth = 0.5, c='#3a5772')          
        ax0.scatter(rs,new_depth,s = s_rad,alpha = 0.7,
                    c='#6291bf',edgecolors='#3a5772')
        
        min = (group.met_cont.min())
        max = (group.met_cont.max()) 
        p = 100*(min/max)
        if p <5:
            p = 0
        p_wat = 100-p
        perc_to_atm.append(p)
        perc_to_wat.append(p_wat)
        mins.append(min)
        maxs.append(max) 
        
    for rad in radii2: 
        group =  gb.get_group(rad)[['met_cont','depth','rad_evol']] 
        new_depth = np.arange(group.depth.min(),group.depth.max(),3)
        f  = interpolate.interp1d(group.depth,group.rad_evol)
        rs = f(new_depth) 
        s_rad = (rs/np.max(rs))*100        
        ax0.plot(group.rad_evol,group.depth,alpha = 0.2,
                 linewidth = 0.5, c='#3a5772')        
        ax0.scatter(rs,new_depth,s = s_rad,alpha = 0.7,
                    c='#e1e1b6',edgecolors='#3a5772')
 
        min = (group.met_cont.min())
        max = (group.met_cont.max()) 
        p = 100*(min/max)
        p_wat = 100-p
        perc_to_atm.append(p)
        perc_to_wat.append(p_wat)
        mins.append(min)
        maxs.append(max)
 
    ax1.bar(radii3,perc_to_atm, width = step,
            bottom=perc_to_wat,edgecolor='k',
            align='center',color = '#e1e1b6',
            label = 'goes to atmosphere')
    
    ax1.bar(radii3,perc_to_wat, width = step,
            align='center', 
            edgecolor='k',color ='#6291bf',
            label = 'dissolves in water')
    
    ax1.set_ylabel('%')
    ax0.set_ylabel('depth, m')
                   
    x_text = - 0.06
    y_text = 1.1    
    ax0.text(x_text, y_text, 'A) ', transform=ax0.transAxes,
                 fontsize=14, fontweight='bold', 
                 va='top', ha='right')    
    ax1.text(x_text, y_text, 'B) ', transform=ax1.transAxes,
                 fontsize=14, fontweight='bold',
                 va='top', ha='right')  
        
    ax0.set_xlabel('bubble radius, mm')    
    ax1.set_xlabel('bubble initial radius, mm')
    plt.legend(facecolor = 'w',frameon= True)
    
    ax0.axis((0,9,80,0))      
    ax1.set_xlim(0,9)
    x_ticks = np.arange(0,9,0.5)   
     
    ax0.set_xticks(x_ticks)
    ax1.set_xticks(x_ticks)

    if save == True:
        plt.savefig('Data/Figure8_Bubbles_Diameters.png')
    else:
        plt.show()


if __name__ == '__main__':     
    plot_rad_evol(True)