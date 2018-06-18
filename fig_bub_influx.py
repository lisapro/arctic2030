import pandas as pd
import matplotlib.pyplot as plt
from glob import glob 
import numpy as np 
import matplotlib.gridspec as gridspec
from scipy import interpolate 
import xarray as xr
from statsmodels.nonparametric.smoothers_lowess import lowess
import statsmodels
import methane_equilibrium as me 

def slblt():
    ds = xr.open_dataset(r"C:\Users\elp\OneDrive\Python_workspace\arctic2030\Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc")          
    d_roms = ds.depth.values 
    temp = ds.temp[70].values # random day!
    sal = ds.sal[70].values
    fg = 1.8*10**(-6) # 1.87987 ppb to check 
    methane2 = []
    for n in range(0,len(temp)):
        t = temp[n]  
        s = sal[n]
        d = d_roms[n]
        met2 = me.calc_methane_depth(t,s,fg,d)[0] /1000000 # to mmol/l
        methane2.append(met2)
    return methane2

##### to read all files. 
def read_all_files():
    files  = glob(r"C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re\total_2160_per1sec_woa_sel_13_decav_*.dat")
    #df_tot = []
    for num,f in enumerate(files): 
        if num+1 == 1:
            df_tot = pd.read_csv(f, delimiter = '  ', header = None,
             names = ['depth',num+1],index_col = 'depth', engine='python') 
            df_tot[num+1] = df_tot[num+1] *1000000
        else: 
            f = pd.read_csv(f, delimiter = '  ', header = None,
            names = ['depth',num+1],index_col = 'depth', engine='python')
            f[num+1] = f[num+1] *1000000
            df_tot[num+1] = f[num+1] 
    df_tot.plot() 
    plt.show()
    
    
def plt_flux(rad):    
    plt.style.use('ggplot')
    path = r"C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re(1)\all_for_sbm_79_mm_tab.dat"
    df = pd.read_csv(path,delimiter = '\t' )
    df4 = df[df.radius == rad].reset_index()
    
     # 3.6 bubble, to milliM
    df4.met_cont = df4.met_cont*1000*3.6
    # 3.6 bubble, to second, to microM
    df4.met_flow = df4.met_flow*1000000*3.6*10 
 
    fig = plt.figure()
    gs = gridspec.GridSpec(1,3)
    ax = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    
    ax.set_title('CH$_4$ dissolution \n$\mu M \cdot sec ^{-1}$')
    ax.plot(df4.met_flow,df4.depth,'--',label = 'r = {}'.format(rad) )
    
    s_flow_mean = np.around(df4.met_flow.mean(),decimals = 2)    
    ax.text(0.3,0.9,'mean ={}'.format(s_flow_mean),transform=ax.transAxes)
    
    ma = np.around(df4.met_cont.max(),decimals = 2)
    mi = np.around(df4.met_cont.min(),decimals = 2)
    to_atm  = np.around(ma - mi,decimals = 2)      
    perc = np.around(mi*100/ma,decimals = 2)  

    ax1.set_title('CH$_4$ in bubbles \n$mM$')    
    ax1.plot(df4.met_cont,df4.depth )  
    ax1.text(0.15,0.9,'To atm = {}\nDissolved={}%'.format(
              to_atm,perc),transform=ax1.transAxes) 
    
    f  = interpolate.interp1d(df4.depth,df4.rad_evol) 
    new_depth = np.arange(1.05,78,4)
    radii = f(new_depth) 
    s_rad = (radii/max(radii))*100 
    

    ax2.set_title('Bubble radius \n$mm$')
    ax2.plot(df4.rad_evol,df4.depth,'-',alpha = 0.1 ) 
    ax2.scatter(radii,new_depth,s = s_rad,alpha = 0.7)

    for axis in (ax,ax1,ax2):
        axis.set_ylim(80,0)
        #axis.legend()
    plt.savefig('Data/bubbles_r_{}.png'.format(rad))    
    #plt.show()
    
#radius(mm) bubble volume(cm^3) methane contnent in bubble (mol) the value distribution density
def init_conditions(): 
    path = r"C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\start_dat_2160 bub_descript.dat"
    df = pd.read_csv(path,delimiter = ' ',header = 0,names = ['radius(mm)','bub_v(cm^3)','met_cont(mol)','dens_dist'],usecols = ['radius(mm)','bub_v(cm^3)','met_cont(mol)','dens_dist']) #,,  engine='python')
    df.plot(subplots = True,kind = 'bar')
    s = df['met_cont(mol)'].sum()
    print (s*1000)
    
def plot_scenario_bbbbb(df4,df3,df2,df1,df_sum,smoothed_flow,new_depth,int_cont):   
     
    # Scenario B1 1 bubble 4mm, 1 bubble 3mm, 1 bubble 2mm and 1 bubble 1mm
    fig = plt.figure()
    gs = gridspec.GridSpec(1,3)
    gs.update(top = 0.8)
    ax = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])    
 
    fig.suptitle('Scenario B1',size = 15)
    # Calculate means,fluxes etc.
    s_flow_mean = np.around(df_sum.met_flow.mean(),decimals = 2)         
    ma = np.around(df_sum.met_cont.max(),decimals = 2)
    mi = np.around(df_sum.met_cont.min(),decimals = 2)
    to_atm  = np.around(ma - mi,decimals = 2)      
    perc = np.around((to_atm*100)/ma,decimals = 2) 
    
    '''labels = ['Dissolved', 'To atmosphere']
    sizes = [perc, 100-perc]
    patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)   ''' 
    
    print ('Mean sum flow from the seep',s_flow_mean)
    print ('To atmosphere ', to_atm)
    print ('Percentage of dissolved CH4 ', perc) 
    
    ax.set_title('CH$_4$ dissolution \n$mM \cdot sec ^{-1}$')
    ax1.set_title('CH$_4$ in bubbles \n$mM$')    
    ax2.set_title('Bubble radius \n$mm$')
        
    l = range(1,5)
    for n,d in enumerate((df1,df2,df3,df4)):
        ax.plot(d.met_flow,d.depth,'--',label = '{}mm'.format(l[n]))
        ax1.plot(d.met_cont,d.depth,'--')   
        ax2.plot(d.rad_evol,d.depth,'-',alpha = 1) 
        


    #ax.plot(int_2_flow,new_depth2,'k--')
    ax.plot(smoothed_flow,new_depth,'k-',label = 'sum')
    #ax1.plot(df_sum.met_cont,df_sum.depth,'k-', label = 'sum') 
    ax1.plot(int_cont,new_depth,'k-', label = 'sum') 
    
    for axis in (ax,ax1,ax2):
        axis.set_ylim(80,0)
    ax.legend()    
    #plt.savefig('Data/Scenario_B1.png')    
    plt.show()
    
def plot_scenario(df_list, df_sum,smoothed_flow,new_depth,int_cont,title):   
     
    # Scenario B2 1 bubble 2mm and 3 bubble 1mm
    
    fig = plt.figure()
    gs = gridspec.GridSpec(1,3)
    gs.update(top = 0.8)
    ax = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])    
 
    fig.suptitle(title,size = 15)
    # Calculate means,fluxes etc.
    s_flow_mean = np.around(df_sum.met_flow.mean(),decimals = 2)         
    ma = np.around(df_sum.met_cont.max(),decimals = 2)
    mi = np.around(df_sum.met_cont.min(),decimals = 2)
    to_atm  = np.around(ma - mi,decimals = 2)  
    print ('max,min', ma,mi)    
    perc = np.around((to_atm*100)/ma,decimals = 2) 
    
    '''labels = ['Dissolved', 'To atmosphere']
    sizes = [perc, 100-perc]
    patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)   ''' 
    
    print ('Mean sum flow from the seep',s_flow_mean)
    print ('To atmosphere ', to_atm)
    print ('Percentage of dissolved CH4 ', perc) 
    
    ax.set_title('CH$_4$ dissolution \n$mM \cdot sec ^{-1}$')
    
    ax1.set_title('CH$_4$ in bubbles \n$mM$')    
    ax2.set_title('Bubble radius \n$mm$')
        
    l = range(1,5)
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
def calculate_scenario_B1(d_roms,pl,days):
    path= r"C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re(1)\all_for_sbm_79_mm_tab.dat" 
    df = pd.read_csv(path,delimiter = '\t' )
    df4 = df[df.radius == 4].reset_index()
    df3 = df[df.radius == 3].reset_index()    
    df2 = df[df.radius == 2].reset_index()      
    df1 = df[df.radius == 1].reset_index()  

    df_sum = df4.loc[:,['depth','rad_evol','met_cont','met_flow','vbub']]

    for col in df_sum.columns[1:] :
        df_sum[col] = df4[col].add(df3[col]).add(
             df2[col]).add(df1[col],fill_value= 0)   
    
    for d in (df1,df2,df3,df4,df_sum):           
        d.met_cont = d.met_cont*1000 # milliM
        d.met_flow = d.met_flow*1000  # milliM/sec
            
    # Calculate means,fluxes etc.
    s_flow_mean = np.around(df_sum.met_flow.mean(),decimals = 2)         
    ma = np.around(df_sum.met_cont.max(),decimals = 2)
    mi = np.around(df_sum.met_cont.min(),decimals = 2)
    to_atm  = np.around(ma - mi,decimals = 2)      
    perc = np.around((to_atm*100)/ma,decimals = 2) 
    print (ma,mi)
    print ('Mean sum flow from the seep',s_flow_mean)
    print ('To atmosphere ', to_atm)
    print ('Percentage of dissolved CH4 ', perc) 

    ## Interpolate to roms depths and 365 days
    new_depth =  d_roms
    f_flow  = interpolate.interp1d(df_sum.depth,
                                   df_sum.met_flow,
                                   fill_value = 'extrapolate',
                                   kind='nearest')  
    f_cont  = interpolate.interp1d(df_sum.depth,
                                   df_sum.met_cont,
                                   fill_value = 'extrapolate',
                                   kind='nearest') 
    int_flow = f_flow(new_depth)
    int_cont = f_cont(new_depth)
    
    smoothed_flow = lowess(new_depth,int_flow,frac=1,is_sorted = False)[:,0]   
    # Reverse array, since in roms axis direction is up    
    cont = np.flip(smoothed_flow,0)
    smoothed_flow = np.flip(smoothed_flow,0)        
    df_flow = pd.DataFrame(index = smoothed_flow,columns = days)
    df_cont = pd.DataFrame(index = cont,columns = days)
    if pl == True: 
        plot_scenario(df4,df3,df2,df1,df_sum,smoothed_flow,new_depth,int_cont)        
        
    for n in days:
        df_flow[n] = smoothed_flow      
        df_cont[n] = cont

    return df_flow.T, df_cont.T


def calculate_baseline(d_roms,days):
    ds = xr.open_dataset(r"C:\Users\elp\OneDrive\Python_workspace\arctic2030\Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc") 
    new_depth =  d_roms
    #days = np.arange(1,367)
    slb = slblt()
    df_slb = pd.DataFrame(index = slb,columns = days)        
    for n in days:
        df_slb[n] = slb     
    return df_slb.T

def calculate_scenario_B1_stop(d_roms,pl,slb,days):
    path= r"C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re(1)\all_for_sbm_79_mm_tab.dat" 
    df = pd.read_csv(path,delimiter = '\t' )
    df4 = df[df.radius == 4].reset_index()
    df3 = df[df.radius == 3].reset_index()    
    df2 = df[df.radius == 2].reset_index()      
    df1 = df[df.radius == 1].reset_index()  

    df_sum = df4.loc[:,['depth','rad_evol','met_cont','met_flow','vbub']]

    for col in df_sum.columns[1:] :
        df_sum[col] = df4[col].add(df3[col]).add(
             df2[col]).add(df1[col],fill_value= 0)   
    
    for d in (df1,df2,df3,df4,df_sum):           
        d.met_cont = d.met_cont*1000 # milliM
        d.met_flow = d.met_flow*1000  # milliM/sec
            
    # Calculate means,fluxes etc.
    s_flow_mean = np.around(df_sum.met_flow.mean(),decimals = 2)         
    ma = np.around(df_sum.met_cont.max(),decimals = 2)
    mi = np.around(df_sum.met_cont.min(),decimals = 2)
    to_atm  = np.around(ma - mi,decimals = 2)      
    perc = np.around((to_atm*100)/ma,decimals = 2) 
    print (ma,mi)
    print ('Mean sum flow from the seep',s_flow_mean)
    print ('To atmosphere ', to_atm)
    print ('Percentage of dissolved CH4 ', perc) 

    ## Interpolate to roms depths and 365 days
    new_depth =  d_roms
    #days = np.arange(1,367)
    f_flow  = interpolate.interp1d(df_sum.depth,
                                   df_sum.met_flow,
                                   fill_value = 'extrapolate',
                                   kind='nearest')  
    f_cont  = interpolate.interp1d(df_sum.depth,
                                   df_sum.met_cont,
                                   fill_value = 'extrapolate',
                                   kind='nearest') 
    int_flow = f_flow(new_depth)
    int_cont = f_cont(new_depth)
    
    smoothed_flow = lowess(new_depth,int_flow,frac=1,is_sorted = False)[:,0]   
    # Reverse array, since in roms axis direction is up    
    cont = np.flip(smoothed_flow,0)
    smoothed_flow = np.flip(smoothed_flow,0)        
    df_flow = pd.DataFrame(index = smoothed_flow,columns = days)
    df_cont = pd.DataFrame(index = cont,columns = days)
    if pl == True: 
        plot_scenario(df4,df3,df2,df1,df_sum,smoothed_flow,new_depth,int_cont)        
        
    for n in days:
        df_flow[n] = slb     
        df_cont[n] = np.nan
    for n in np.arange(1,150):
        df_flow[n] = smoothed_flow      
        df_cont[n] = cont
    return df_flow.T, df_cont.T




def calculate_scenario_B1_50(d_roms,pl,slb,days):
    ## Sampe seep as B1 but working 50% of time 
    path= r"C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re(1)\all_for_sbm_79_mm_tab.dat" 
    df = pd.read_csv(path,delimiter = '\t' )
    df4 = df[df.radius == 4].reset_index()
    df3 = df[df.radius == 3].reset_index()    
    df2 = df[df.radius == 2].reset_index()      
    df1 = df[df.radius == 1].reset_index()  

    df_sum = df4.loc[:,['depth','rad_evol','met_cont','met_flow','vbub']]

    for col in df_sum.columns[1:] :
        df_sum[col] = df4[col].add(df3[col]).add(
             df2[col]).add(df1[col],fill_value= 0)   
    
    for d in (df1,df2,df3,df4,df_sum):           
        d.met_cont = d.met_cont*1000 # milliM
        d.met_flow = d.met_flow*1000  # milliM/sec
            
    # Calculate means,fluxes etc.
    s_flow_mean = np.around(df_sum.met_flow.mean(),decimals = 2)         
    ma = np.around(df_sum.met_cont.max(),decimals = 2)
    mi = np.around(df_sum.met_cont.min(),decimals = 2)
    to_atm  = np.around(ma - mi,decimals = 2)      
    perc = np.around(mi*100/ma,decimals = 2) 
    
    print ('Mean sum flow from the seep',s_flow_mean)
    print ('To atmosphere ', to_atm)
    print ('Percentage of dissolved CH4 ', perc) 

    ## Interpolate to roms depths and 365 days
    new_depth =  d_roms
    #days = np.arange(1,367)
    f_flow  = interpolate.interp1d(df_sum.depth,
                                   df_sum.met_flow,
                                   fill_value = 'extrapolate',
                                   kind='nearest')  
    f_cont  = interpolate.interp1d(df_sum.depth,
                                   df_sum.met_cont,
                                   fill_value = 'extrapolate',
                                   kind='nearest') 
    int_flow = f_flow(new_depth)
    int_cont = f_cont(new_depth)
    
    smoothed_flow = lowess(new_depth,int_flow,frac=1,is_sorted = False)[:,0]   
    # Reverse array, since in roms axis direction is up    
    cont = np.flip(smoothed_flow,0)
    smoothed_flow = np.flip(smoothed_flow,0)        
    df_flow = pd.DataFrame(index = smoothed_flow,columns = days)
    df_cont = pd.DataFrame(index = cont,columns = days)
    if pl == True: 
        df_list = [df4,df3,df2,df1]
        plot_scenario(df_list,df_sum,smoothed_flow,new_depth,int_cont,'Scenario B1 50%')        
    
        
    for n in days:
        df_flow[n] = slb #smoothed_flow      
        df_cont[n] = np.nan #cont

    import random
    lst = range(1,366)
    n = int(366/2)
    rand = random.sample(lst,n)

    for n in rand:
        df_flow[n] = smoothed_flow      
        df_cont[n] = cont
                
    return df_flow.T, df_cont.T
#init_conditions() 

def calculate_scenario_B2(d_roms,pl,days):
    path= r"C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re(1)\all_for_sbm_79_mm_tab.dat" 
    df = pd.read_csv(path,delimiter = '\t' )
    
    #df3 = df[df.radius == 3].reset_index()    
    df2 = df[df.radius == 2].reset_index()          
    df1 = df[df.radius == 1].reset_index()  
    
    df_sum = df2.loc[:,['depth','rad_evol','met_cont','met_flow','vbub']]
    
    # Tree bubbles with d 1 (so we add df1 three times to df_sum) 
    for col in df_sum.columns[1:] :
        df_sum[col] = df2[col].add(df1[col],fill_value= 0).add(df1[col],fill_value= 0).add(df1[col],fill_value= 0)     
    
    for d in (df1,df2,df_sum):           
        d.met_cont = d.met_cont*1000000 # milliM
        d.met_flow = d.met_flow*1000000  # milliM/sec
            
    # Calculate means,fluxes etc.
    s_flow_mean = np.around(df_sum.met_flow.mean(),decimals = 2)         
    ma = np.around(df_sum.met_cont.max(),decimals = 2)
    mi = np.around(df_sum.met_cont.min(),decimals = 2)
    to_atm  = np.around(ma - mi,decimals = 2)      
    perc = np.around(mi*100/ma,decimals = 2) 
    
    print ('Mean sum flow from the seep',s_flow_mean)
    print ('To atmosphere ', to_atm)
    print ('Percentage of dissolved CH4 ', perc) 

    ## Interpolate to roms depths and 365 days
    new_depth =  d_roms
    #days = np.arange(1,367)
    f_flow  = interpolate.interp1d(df_sum.depth,
                                   df_sum.met_flow,
                                   fill_value = 'extrapolate',
                                   kind='nearest')  
    f_cont  = interpolate.interp1d(df_sum.depth,
                                   df_sum.met_cont,
                                   fill_value = 'extrapolate',
                                   kind='nearest') 
    int_flow = f_flow(new_depth)
    int_cont = f_cont(new_depth)
    
    smoothed_flow = lowess(new_depth,int_flow,frac=1,is_sorted = False)[:,0]   
    # Reverse array, since in roms axis direction is up    
    cont = np.flip(smoothed_flow,0)
    smoothed_flow = np.flip(smoothed_flow,0)        
    df_flow = pd.DataFrame(index = smoothed_flow,columns = days)
    df_cont = pd.DataFrame(index = cont,columns = days)
    if pl == True: 
        df_list =  [df2,df1]
        plot_scenario(df_list,df_sum,smoothed_flow,new_depth,int_cont,'Scenario B2')        
        
    for n in days:
        df_flow[n] = smoothed_flow      
        df_cont[n] = cont

    return df_flow.T, df_cont.T

def calculate_scenario_B_10min(d_roms,pl):
    '''mu, sigma = 4, 0.5 # mean and standard deviation
    sizes =    sorted(np.round(np.random.normal(mu, sigma, 2160),decimals = 2))     
    print (sizes)
    plt.hist(sizes)'''
    #count, bins, ignored = plt.hist(sizes, 30, normed=True)
    #plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
    #            np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
    #      linewidth=2, color='r')
    plt.show()

    path= r"C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re(1)\all_for_sbm_79_mm_tab.dat" 
    df = pd.read_csv(path,delimiter = '\t' )
    
    #df3 = df[df.radius == 3].reset_index()    
    df2 = df[df.radius == 2].reset_index()          
    df1 = df[df.radius == 1].reset_index()  
    df_all = pd.merge(df1,df2,on=['depth','time'],suffixes = ['_1mm','_2mm'])
    print (df_all.head())
    df_sum = df2.loc[:,['depth','rad_evol','met_cont','met_flow','vbub']]
    
    # Tree bubbles with d 1 (so we add df1 three times to df_sum) 
    for col in df_sum.columns[1:] :
        df_sum[col] = df2[col].add(df1[col],fill_value= 0).add(df1[col],fill_value= 0).add(df1[col],fill_value= 0)     
    
    for d in (df1,df2,df_sum):           
        d.met_cont = d.met_cont*1000000 # milliM
        d.met_flow = d.met_flow*1000000  # milliM/sec
            
    # Calculate means,fluxes etc.
    s_flow_mean = np.around(df_sum.met_flow.mean(),decimals = 2)         
    ma = np.around(df_sum.met_cont.max(),decimals = 2)
    mi = np.around(df_sum.met_cont.min(),decimals = 2)
    to_atm  = np.around(ma - mi,decimals = 2)      
    perc = np.around(mi*100/ma,decimals = 2) 
    
    print ('Mean sum flow from the seep',s_flow_mean)
    print ('To atmosphere ', to_atm)
    print ('Percentage of dissolved CH4 ', perc) 

    ## Interpolate to roms depths and 365 days
    new_depth =  d_roms
    days = np.arange(1,367)
    f_flow  = interpolate.interp1d(df_sum.depth,
                                   df_sum.met_flow,
                                   fill_value = 'extrapolate',
                                   kind='nearest')  
    f_cont  = interpolate.interp1d(df_sum.depth,
                                   df_sum.met_cont,
                                   fill_value = 'extrapolate',
                                   kind='nearest') 
    int_flow = f_flow(new_depth)
    int_cont = f_cont(new_depth)
    
    smoothed_flow = lowess(new_depth,int_flow,frac=1,is_sorted = False)[:,0]   
    # Reverse array, since in roms axis direction is up    
    cont = np.flip(smoothed_flow,0)
    smoothed_flow = np.flip(smoothed_flow,0)        
    df_flow = pd.DataFrame(index = smoothed_flow,columns = days)
    df_cont = pd.DataFrame(index = cont,columns = days)
    if pl == True: 
        df_list =  [df2,df1]
        plot_scenario(df_list,df_sum,smoothed_flow,new_depth,int_cont,'Scenario B2')        
        
    for n in days:
        df_flow[n] = smoothed_flow      
        df_cont[n] = cont

    return df_flow.T, df_cont.T 


if __name__ == '__main__':     
    ds = xr.open_dataset(r"C:\Users\elp\OneDrive\Python_workspace\arctic2030\Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc")          
    roms_levels = ds.depth.values 

    #calculate_scenario_B1_50(roms_levels,True,solubility)
    calculate_baseline(ds,roms_levels)
    #calculate_scenario_B2(roms_levels,True)
    #calculate_scenario_B_10min(roms_levels,True)
    
    #plot_scenario_B1() 
    #read_all_files()
    ##plot_scenario_B1(roms_levels) 