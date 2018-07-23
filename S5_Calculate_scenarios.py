import pandas as pd 
import numpy as np       
from scipy import interpolate 
from statsmodels.nonparametric.smoothers_lowess import lowess

def make_df_sum(s,scen):
    path= r"C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re(1)\all_for_sbm_79_mm_tab.dat" 
    df = pd.read_csv(path,delimiter = '\t' )    
    sizes,fraction = s
    df1 = df[df.radius == sizes[0]].reset_index()
    df_sum = df1.loc[:,['depth','rad_evol','met_cont','met_flow','vbub']]
    
    for n,num in enumerate(sizes):          
        df2 = df[df.radius == n].reset_index()      
        for col in df_sum.columns[1:] :
            df_sum[col] = df_sum[col].add(df2[col],fill_value= 0)   
              
    df_sum.met_cont = df_sum.met_cont*1000 # milliM
    df_sum.met_flow = df_sum.met_flow*1000  # milliM/m2/sec  
    if scen == 'B2_10_30min':
          df_sum.met_flow = df_sum.met_flow*(30/24*60)
    return df_sum

def calculate_scenarios(d_roms,pl,days,sc): #slb,
    # get sum flux for the scenario
    # bubble sizes, fraction 3 is 1/3 
    scenario_dict = {'B0_30':[[4,2,1],0.3],
                     'B2_30':[[2,1],0.3],
                     'B2_10':[[2,1],0.1],
                     'B2_10_30min':[[2,1],0.1],
                     'B0_50':[[4,2,1],0.5],
                     'B1_30':[[4,3,2,1],0.3],
                     'B1_50':[[4,3,2,1],0.5],
                     'B1':[[4,3,2,1],1]}
    
    df_sum = make_df_sum(scenario_dict[sc],sc) 
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
    print ("new depth",new_depth)
    int_flow = f_flow(new_depth)
    int_cont = f_cont(new_depth)
    
    smoothed_flow = lowess(new_depth,int_flow,frac=1,is_sorted = False)[:,0]   
    
    # Reverse array, since in roms axis direction is up    
    cont = np.flip(smoothed_flow,0)
    smoothed_flow = np.flip(smoothed_flow,0)      
    smoothed_flow_ice = smoothed_flow.copy()
    smoothed_flow_ice[-1:] = to_atm  

    df_flow = pd.DataFrame(index = smoothed_flow,columns = days)
    df_cont = pd.DataFrame(index = cont,columns = days)
    
    if frac ==1 :
        for n in days:
            if n < 215 or n> 308: # All the rest methane during ice-covered season
                df_flow[n] = smoothed_flow_ice
                df_cont[n] = cont
            else:
                df_flow[n] = smoothed_flow      
                df_cont[n] = cont               
    
    else:        
        for n in days:
            df_flow[n] = 0 #slb #smoothed_flow      
            df_cont[n] = np.nan #cont

        import random
        lst = range(1,366)
        n = int(365*frac)
        rand = random.sample(lst,n)
    
        for n in rand:
            if n < 215 or n> 308: 
                df_flow[n] = smoothed_flow_ice
                df_cont[n] = cont
            else:
                df_flow[n] = smoothed_flow      
                df_cont[n] = cont    

    return df_flow.T, df_cont.T