import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np 
from scipy import interpolate 
import pandas as pd


def plt_flux_1bub(rad,save = False):    
    '''
    Function plots the dissolution rate ,
    methane content in the bubble, and 
    bubble radius evolutions 
    based on calculations from the
    Single bubble model 
    https://github.com/zagoven/bubl
    '''
    
    plt.style.use('ggplot')
    path = r"C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re(1)\all_for_sbm_79_mm_tab.dat"
    df = pd.read_csv(path,delimiter = '\t' )
    df2 = df[df.radius == rad].reset_index()

    df2.met_cont = df2.met_cont*1000  
    df2.met_flow = df2.met_flow*1000000   # to microM/second 
    
    fig = plt.figure(figsize =(6,4), dpi=100)
    gs = gridspec.GridSpec(1,3)
    gs.update(top = 0.8)
    ax = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    ax1.set_yticklabels([])
    ax2.set_yticklabels([])
    ax.set_title('CH$_4$ dissolution \n$\mu M \cdot sec ^{-1}$')
    ax.plot(df2.met_flow,df2.depth,'--',label = 'r = {}'.format(2))
    
    s_flow_mean = np.around(df2.met_flow.mean(),decimals = 2)   
        
    ma = np.around(df2.met_cont.max(),decimals = 4)
    mi = np.around(df2.met_cont.min(),decimals = 4)    
    to_atm  = np.around(ma - mi,decimals = 4)      
    perc = np.around(mi*100/ma,decimals = 4)  
    
    print (ma,mi,to_atm,perc)
       
    ax1.set_title('CH$_4$ in bubbles \n$mM$')    
    ax1.plot(df2.met_cont,df2.depth )  
     
    f  = interpolate.interp1d(df2.depth,df2.rad_evol)

    new_depth = np.arange(1.05,78,4)
    radii = f(new_depth) 
    s_rad = (radii/max(radii))*100 

    ax2.set_title('Bubble radius \n$mm$')
    ax2.plot(df2.rad_evol,df2.depth,'-',alpha = 0.1 )  
    ax2.scatter(radii,new_depth,s = s_rad,alpha = 0.7)
    
    for axis in (ax,ax1,ax2):
        axis.set_ylim(80,0)        
        
    if save == True:
        plt.savefig('Data/bubbles_r_{}.png'.format(rad),transparent = False) 
    else:  
        plt.show()    
        
if __name__ == '__main__':   
    plt_flux_1bub(4) 
 