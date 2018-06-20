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

def calculate_baseline(d_roms,days):
    ds = xr.open_dataset(r"C:\Users\elp\OneDrive\Python_workspace\arctic2030\Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc") 
    new_depth =  d_roms
    #days = np.arange(1,367)
    slb = slblt()
    df_slb = pd.DataFrame(index = slb,columns = days)        
    for n in days:
        df_slb[n] = slb     
    return df_slb.T


def make_df_sum(s):
    path= r"C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re(1)\all_for_sbm_79_mm_tab.dat" 
    df = pd.read_csv(path,delimiter = '\t' )    
    sizes,fraction = s
    print(sizes)
    df1 = df[df.radius == sizes[0]].reset_index()
    df_sum = df1.loc[:,['depth','rad_evol','met_cont','met_flow','vbub']]
    
    for n,num in enumerate(sizes):          
        df2 = df[df.radius == n].reset_index()      
        for col in df_sum.columns[1:] :
            df_sum[col] = df_sum[col].add(df2[col],fill_value= 0)   
              
    df_sum.met_cont = df_sum.met_cont*1000 # milliM
    df_sum.met_flow = df_sum.met_flow*1000  # milliM/sec    
    return df_sum


def calculate_scenarios(d_roms,pl,slb,days,sc):
    # get sum flux for the scenario
    # bubble sizes, fraction 3 is 1/3 
    scenario_dict = {'B0_30':[[4,2,1],3],
                     'B2_30':[[2,1],3],
                     'B0_50':[[4,2,1],2],
                     'B1_30':[[4,3,2,1],3],
                     'B1_50':[[4,3,2,1],2],
                     'B1':[[4,3,2,1],1]}
    df_sum = make_df_sum(scenario_dict[sc]) 
    frac = scenario_dict[sc][1]       
    # Calculate means,fluxes etc.
    s_flow_mean = np.around(df_sum.met_flow.mean(),decimals = 2)         
    ma = np.around(df_sum.met_cont.max(),decimals = 2)
    mi = np.around(df_sum.met_cont.min(),decimals = 2)
    to_atm  = np.around(ma - mi,decimals = 2)      
    perc = np.around(mi*100/ma,decimals = 2) 
    perc1 = np.around(100-perc,decimals = 2)
    print ('Mean sum flow from the seep {}'.format(sc),s_flow_mean)
    print ('To atmosphere  {} '.format(sc), to_atm)
    print ('Percentage of dissolved CH4 {}'.format(sc), perc) 
    print ('Percentage of flux CH4 to atm {}'.format(sc), perc1) 
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
    
    '''plt.plot(int_flow,new_depth,'--')
    plt.plot(smoothed_flow,new_depth)
    plt.ylim(80,0)
    plt.show()  '''   
    if frac ==1 :
        for n in days:
            df_flow[n] = smoothed_flow      
            df_cont[n] = cont        
    else:           
        for n in days:
            df_flow[n] = slb #smoothed_flow      
            df_cont[n] = np.nan #cont
        
        import random
        lst = range(1,366)
        n = int(365/frac)
        rand = random.sample(lst,n)
    
        for n in rand:
            df_flow[n] = smoothed_flow      
            df_cont[n] = cont

    return df_flow.T, df_cont.T




if __name__ == '__main__':     
    pass
    #ds = xr.open_dataset(r"C:\Users\elp\OneDrive\Python_workspace\arctic2030\Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc")          
    #roms_levels = ds.depth.values 
    #days = np.arange(1,366)
    #slb_year = calculate_baseline(roms_levels,days)
    #df_flow = calculate_scenarios(roms_levels,True,slb_year,days,'B0_30')[0]

    #calculate_scenario_B1_50(roms_levels,True,solubility)
    #calculate_baseline(ds,roms_levels)
    #calculate_scenario_B2(roms_levels,True)
    #calculate_scenario_B_10min(roms_levels,True)
    
    #plot_scenario_B1() 
    #read_all_files()
    ##plot_scenario_B1(roms_levels) 