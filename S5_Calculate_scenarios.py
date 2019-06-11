import pandas as pd 
import numpy as np       
from scipy import interpolate 
from statsmodels.nonparametric.smoothers_lowess import lowess
import xarray as xr

global roms_path,bub_path
roms_path = r"Data/ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc"
bub_path = r'Data/all_for_sbm_79_mm_tab.dat'

def make_df_sum(s,scen):
    ''' Create scenario, sum fluxes from 
        different bubles '''
    df = pd.read_csv(bub_path,delimiter = '\t' )    
    sizes,fraction = s
    df1 = df[df.radius == sizes[0]].reset_index()
    df_sum = df1.loc[:,['depth','rad_evol','met_cont','met_flow','vbub']]
    
    for n,num in enumerate(sizes):          
        df2 = df[df.radius == num].reset_index()      
        for col in df_sum.columns[1:] :
            df_sum[col] = df_sum[col].add(df2[col],fill_value= 0)   
                      
    df_sum.met_cont = df_sum.met_cont*1000 # milliM
    df_sum.met_flow = df_sum.met_flow*1000 # milliM/m2/sec  
    return df_sum

def calculate_spin_up(z,days):
    return pd.DataFrame(0, index = z,columns = days).T

def slblt_1d():
    import S1_methane_equilibrium  as me
    ds = xr.open_dataset(roms_path)      
    d_roms = ds.depth.values 
    temp = ds.temp[70].values # random day!
    sal = ds.sal[70].values
    fg = 1.8*10**(-6) # 1.87987 ppb to check 
    methane2 = []
    for n in range(0,len(temp)):
        met2 = me.calc_methane_depth(temp[n],sal[n],fg,d_roms[n])[0] /1.e6 # to mmol/l
        methane2.append(met2)
    return methane2


def calculate_baseline(days):
    import S1_methane_equilibrium  as me
    df_slb = me.calculate_equilibrium_solubility(days) 
    return df_slb #.T

def calculate_scenarios(d_roms,days,sc): 
    ''' get sum flux for the scenario
       bubble sizes'''
    scenario_dict = {#'B0_30':[[4,2,1],0.3],
                     #'B2_50':[[2,1],0.5],
                     #'B3_50':[[6,3,2],0.5],                     
                     'BOM_Basic_seep':[[4]*4,0.5],                    
                     #'O_Increased_oxidation_rate':[[4]*4,0.5],                    
                     #'M_Increased_mixing_rate':[[4]*4,0.5],
                     'S_Small_bubbles':[[2]*32,0.5],
                     'F_Reduced_flux':[[4]*3,0.5],
                     'F2_Reduced_flux':[[4]*2,0.5]}
    
    df_sum = make_df_sum(s = scenario_dict[sc],scen = sc)
    frac = scenario_dict[sc][1]    

    # Calculate means,fluxes etc.    
    ma = df_sum.met_cont.max()
    mi = df_sum.met_cont.min()
    to_atm  = np.around(mi,decimals=2) 
    perc = np.around((mi/ma)*100,decimals=2)
    perc1 = np.around(100-perc,decimals=2)
    
    print ('Sum flow to water milliM/sec from the whole! seep {}'.format(sc),df_sum.met_flow.sum())    
    #print ('Mean flow to water milliM/m2/sec from one horizont {}'.format(sc),df_sum.met_flow.mean())
    print ('To atmosphere  {} # milliM/sec '.format(sc), to_atm, 'max',ma,'min',mi)
    print ('CH4 dissolved in the water during uplifting,scenario:{} {} %'.format(sc, perc1)) 
    print ('Percentage of flux CH4 to atm {}  {}%'.format(sc, perc)) 
    print ("     ")     

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
    #int_cont = f_cont(new_depth)
    
    smoothed_flow = lowess(new_depth,int_flow,frac=1,is_sorted = False)[:,0]   
    
    # Reverse array, since in roms axis direction is up    
    cont = np.flip(smoothed_flow,0)
    smoothed_flow = np.flip(smoothed_flow,0)      
    smoothed_flow_ice = smoothed_flow.copy()
    smoothed_flow_ice[-1:] = to_atm  

    df_flow = pd.DataFrame(index = new_depth,columns = days)
    df_cont = pd.DataFrame(index = new_depth,columns = days)
    
    start_wat = 236 # 27 aug
    stop_wat = 303 # 28 oct 
    
    if frac ==1 :
        for n in days:
            if n < start_wat or n> stop_wat: 
                df_flow[n] = smoothed_flow_ice
                df_cont[n] = cont
            else:
                df_flow[n] = smoothed_flow  # All the rest methane during ice-covered season    
                df_cont[n] = cont                   
    else:        
        for n in days:
            df_flow[n] = 0 #slb #smoothed_flow      
            df_cont[n] = np.nan #cont
        import random
        random.seed(14)
        lst = range(1,366)
        n = int(365*frac)
        rand = random.sample(lst,n)

        for n in rand:
            if n < start_wat or n> stop_wat: 
                df_flow[n] = smoothed_flow_ice
                df_cont[n] = cont
            else:
                df_flow[n] = smoothed_flow      
                df_cont[n] = cont    

    return df_flow.T #, df_cont.T

if __name__ == '__main__': 
    pass
