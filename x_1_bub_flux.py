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

def plot_scenario(df_list, df_sum,smoothed_flow,
                  new_depth,int_cont,title):   
     
    ''' Scenario B2 1 bubble 2mm and 3 bubble 1mm '''
    fig = plt.figure()
    gs = gridspec.GridSpec(1,3)
    gs.update(top = 0.8)
    ax = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])    
    fig.suptitle(title,size = 15)
    
    s_flow_mean = np.around(df_sum.met_flow.mean(),decimals = 2)         
    ma = np.around(df_sum.met_cont.max(),decimals = 2)
    mi = np.around(df_sum.met_cont.min(),decimals = 2)
    to_atm  = np.around(ma - mi,decimals = 2)  
    perc = np.around((to_atm*100)/ma,decimals = 2) 
    
    print ('max,min', ma,mi)
    print ('Mean sum flow from the seep',s_flow_mean)
    print ('To atmosphere ', to_atm)
    print ('Percentage of dissolved CH4 ', perc) 
    
    ax.set_title('CH$_4$ dissolution \n$mM \cdot sec ^{-1}$')    
    ax1.set_title('CH$_4$ in bubbles \n$mM$')    
    ax2.set_title('Bubble radius \n$mm$')
        
    for n,d in enumerate(df_list):
        ax.plot(d.met_flow,d.depth,'--',label = '{}mm'.format(l[n]))
        ax1.plot(d.met_cont,d.depth,'--')   
        ax2.plot(d.rad_evol,d.depth,'-',alpha = 1) 
        
    #ax.plot(int_2_flow,new_depth2,'k--')
    ax.plot(smoothed_flow,new_depth,'k-',label = 'sum')
    #ax.ticklabel_format(style='sci',scilimits=(0.00001,0.0001),axis='both')
    #ax1.plot(df_sum.met_cont,df_sum.depth,'k-', label = 'sum') 
    ax1.plot(int_cont,new_depth,'k-', label = 'sum') 
    
    for axis in (ax,ax1,ax2):
        axis.set_ylim(80,0)
    ax.legend()    
    #plt.savefig('Data/Scenario_B1.png')    
    plt.show() 

def plt_flux(save = False):   
    ''' 'Figure for 2 bubble/sec scenario  '''
    path = r"C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re(1)\all_for_sbm_79_mm_tab.dat"
    df = pd.read_csv(path,delimiter = '\t' )
    df = df.where(df.met_cont >= 0, 0)

    df2 = df[df.radius == 2].reset_index()
    df1 = df[df.radius == 1].reset_index()
     
    df2.met_cont = df2.met_cont*1000 # milliM
    df1.met_cont = df1.met_cont*1000   
    df2.met_flow = df2.met_flow*1000000 
    df1.met_flow = df1.met_flow*1000000 # microM
    
    fig = plt.figure(figsize =(5,3), dpi=100)
    gs = gridspec.GridSpec(1,3)
    gs.update(top = 0.8)
    ax = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    ax1.set_yticklabels([])
    ax2.set_yticklabels([])
    ax.set_title('CH$_4$ dissolution \n$\mu M \cdot sec ^{-1}$')
    ax.plot(df2.met_flow,df2.depth,'--', label = 'r = {}'.format(2))
    ax.plot(df1.met_flow,df1.depth,'--', label = 'r = {}'.format(1))
    
    s_flow_mean = np.around(df2.met_flow.mean(),decimals = 2)    
    
    ma = np.around(df2.met_cont.max(),decimals = 4)
    mi = np.around(df2.met_cont.min(),decimals = 4)
    
    to_atm  = np.around(ma - mi,decimals = 4)      
    perc = np.around(mi*100/ma,decimals = 4)  
    
    print ('Max', ma,'Min',mi,to_atm,perc)
    
    ax1.set_title('CH$_4$ in bubbles \n$mM$')    
    ax1.plot(df2.met_cont,df2.depth,'--' )  
    ax1.plot(df1.met_cont,df1.depth,'--' )  
    
    ax2.set_title('Bubble radius \n$mm$')
    ax2.plot(df2.rad_evol,df2.depth,'-',alpha = 0.1 ) 
    ax2.plot(df1.rad_evol,df1.depth,'-',alpha = 0.1 )    
    
    f  = interpolate.interp1d(df2.depth,df2.rad_evol)
    f2 = interpolate.interp1d(df1.depth,df1.rad_evol) 
    
    new_depth = np.arange(1.05,79,4)
    new_depth2 = np.arange(58,79,1)
    
    radii = f(new_depth) 
    radii2 = f2(new_depth2) 
    
    s_rad = (radii/max(radii))*100 
    s_rad2 = (radii2/max(radii2))*100 
    
    ax2.scatter(radii,new_depth,s = s_rad,alpha = 0.7)
    ax2.scatter(radii2,new_depth2,s = s_rad2,alpha = 0.7)
    
    [axis.set_ylim(80,0) for axis in (ax,ax1,ax2)]    
       
    if save == True: 
        plt.savefig('Data/bubbles_r_B0.png',transparent = True)    
    else:   
        plt.show()
        
                
if __name__ == '__main__':   
    #plt_flux_1bub(4) 
    plt_flux()