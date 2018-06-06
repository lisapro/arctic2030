import pandas as pd
import matplotlib.pyplot as plt
from glob import glob 
import numpy as np 
import matplotlib.gridspec as gridspec
from scipy import interpolate 
    
##### to read all files. 
def read_all_files():
    files  = glob(r'C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re\total_2160_per1sec_woa_sel_13_decav_*.dat')
    #df_tot = []
    for num,f in enumerate(files): 
        if num+1 == 1:
            df_tot = pd.read_csv(f, delimiter = '  ', header = None, names = ['depth',num+1],index_col = 'depth', engine='python') 
            df_tot[num+1] = df_tot[num+1] *1000000
        else: 
            f = pd.read_csv(f, delimiter = '  ', header = None, names = ['depth',num+1],index_col = 'depth', engine='python')
            f[num+1] = f[num+1] *1000000
            #print (f)
            df_tot[num+1] = f[num+1]
        #df_tot.append(pd.read_csv(f, delimiter = '  ', header = None, names = ['depth',num+1],index_col = 'depth', engine='python')) 
    df_tot.plot() 
    #print (df_tot.head())
    plt.show()
    
    
def plt_flux(rad):    

    #'''file  = (r'C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re\total_2160_per1sec_woa_sel_13_decav_05.dat')
    #df = pd.read_csv(file, delimiter = '  ', header = None, names = ['depth','met'], engine='python')
    #df = df[df.depth<79]
    #df.met = df.met * 1000 #000 # to micromoles
    #plt.plot(df.met,df.depth)
    #plt.ylim(80,0)'''
    plt.style.use('ggplot')
    path = r'C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re (1)\all_for_sbm_79_mm_tab.dat'
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
    path = r'C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\start_dat_2160 bub_descript.dat'
    df = pd.read_csv(path,delimiter = ' ',header = 0,names = ['radius(mm)','bub_v(cm^3)','met_cont(mol)','dens_dist'],usecols = ['radius(mm)','bub_v(cm^3)','met_cont(mol)','dens_dist']) #,,  engine='python')
    df.plot(subplots = True,kind = 'bar')
    #plt.show()
    s = df['met_cont(mol)'].sum()
    print (s*1000)
    #print (df)
    
#plt_flux(3)
    
def scenario_B1():    
    # Scenario B1 1 bubble 4mm, 1 bubble 3mm, 1 bubble 2mm and 1 bubble 1mm
    
    plt.style.use('ggplot')
    path = r'C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re (1)\all_for_sbm_79_mm_tab.dat'
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
        # 1 bubble, to milliM
        d.met_cont = d.met_cont*1000 #*3.6
        # 1 bubble, to second, to microM
        d.met_flow = d.met_flow*1000000*10 

    fig = plt.figure()
    gs = gridspec.GridSpec(1,3)
    ax = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])    
 
 
    # Calculate means,fluxes etc.
    s_flow_mean = np.around(df_sum.met_flow.mean(),decimals = 2)         
    ma = np.around(df_sum.met_cont.max(),decimals = 2)
    mi = np.around(df_sum.met_cont.min(),decimals = 2)
    to_atm  = np.around(ma - mi,decimals = 2)      
    perc = np.around(mi*100/ma,decimals = 2) 
    print ('Mean sum flow from the seep',s_flow_mean)
    print ('To atmosphere ', to_atm)
    print ('Percentage of dissolved CH4 ', perc) 
    
    ax.set_title('CH$_4$ dissolution \n$\mu M \cdot sec ^{-1}$')
    ax1.set_title('CH$_4$ in bubbles \n$mM$')    
    ax2.set_title('Bubble radius \n$mm$')
        
    l = range(1,5)
    for n,d in enumerate((df1,df2,df3,df4)):
        ax.plot(d.met_flow,d.depth,'--',label = '{}mm'.format(l[n]))
        ax1.plot(d.met_cont,d.depth,'--')   
        ax2.plot(d.rad_evol,d.depth,'-',alpha = 1) 
        
           
    #ax.plot(df_sum.met_flow,df_sum.depth,label = 'sum')
    
    new_depth = np.arange(0,78,0.01) 
    new_depth2 = np.arange(1,78,1) 
    f_flow  = interpolate.interp1d(df_sum.depth,df_sum.met_flow, fill_value = 'extrapolate',kind='nearest')  
    int_flow = f_flow(new_depth)
    
    #s = interpolate.splrep(int_flow, new_depth, s=0) 
    #int2_flow = s(new_depth)
    #print (s)
    #ax.plot(int_2_flow,new_depth2,'k--')
    #ax.plot(int2_flow,new_depth,'r--')
    
    ax1.plot(df_sum.met_cont,df_sum.depth) 
    
    
                  
    #ax.text(0.3,0.9,'mean ={}'.format(s_flow_mean),transform=ax.transAxes)                  
    #ax1.text(0.15,0.9,'To atm = {}\nDissolved={}%'.format(
    #          to_atm,perc),transform=ax1.transAxes) 
    
    '''
    f  = interpolate.interp1d(df_sum.depth,df_sum.rad_evol) 
    
    radii = f(new_depth) 
    s_rad = (radii/max(radii))*100     
    ax2.scatter(radii,new_depth,s = s_rad,alpha = 0.7)
    '''
    for axis in (ax,ax1,ax2):
        axis.set_ylim(80,0)
    ax.legend()
    
    #plt.savefig('Data/Scanerio_B1_r_{}.png'.format(rad))  
      
    plt.show()
#init_conditions() 
scenario_B1()      
read_all_files()