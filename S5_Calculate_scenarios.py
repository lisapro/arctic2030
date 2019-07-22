import pandas as pd 
import numpy as np       
from scipy import interpolate 
from statsmodels.nonparametric.smoothers_lowess import lowess
import xarray as xr

global roms_path,bub_path
roms_path = r"Data/ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc"
bub_path = r'Data/all_for_sbm_79_mm_tab.dat'

def make_df_sum(s):
    ''' Create scenario, sum rates of dissolution  from 
        different bubbles milliM/sec '''

    # Get rates of dissolution calculated in the Single Bubble Model
    df = pd.read_csv(bub_path,delimiter = '\t' )    

    sizes = s[0]
    df1 = df[df.radius == sizes[0]].reset_index()

    df_sum = df1.loc[:,['depth','rad_evol','met_cont','met_flow','vbub']]
    

    # Make sum according to scenario 
    for n,num in enumerate(sizes):          
        df2 = df[df.radius == num].reset_index()      
        for col in df_sum.columns[1:] :
            df_sum[col] = df_sum[col].add(df2[col],fill_value= 0)   

    #frac = s[1]  # fraction of the day the seep is active 
    df_sum.met_cont = df_sum.met_cont*1000 #* frac # milliM
    df_sum.met_flow = df_sum.met_flow*1000 #* frac  # milliM/sec  #50% of the time 
    return df_sum #Rates of dissolution for sum of bubbles for each scenario, milliM/sec

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

def rate_to_flow(smoothed_flow,new_depth):
    # smoothed_flow # milliM/sec  divide by volume of cell #
    difs_x = [(np.abs(new_depth[n] - new_depth[n-1]) + np.abs(new_depth[n+1] - new_depth[n]))/2. for n in range(1,len(new_depth)-1)]
    last = (np.abs(new_depth[-2] - new_depth[-1]) +  new_depth[-1])/2 #top layer
    first = np.abs(new_depth[1] - new_depth[0])    
    difs_x = np.append(difs_x,last)
    difs_x = np.append(first,difs_x)
    volumes =  10 * difs_x  # volume of cells m3
    new_flow = np.array(smoothed_flow) / np.array(volumes)
    return new_flow,difs_x,volumes

def word_table(table,sc):
    import docx   
    doc = docx.Document()
    doc_table = doc.add_table(rows=(table.shape[0]+1), cols=table.shape[1])  # First row are table headers
    doc_table.style = 'Table Grid'

    for j in range(table.shape[-1]):
        doc_table.cell(0,j).text = table.columns[j]

    for i in range(table.shape[0]):
        for j in range(table.shape[-1]):
            doc_table.cell(i+1,j).text = str("{:.2E}".format(table.values[i,j]))

    return doc.save('./dataframe_{}.docx'.format(sc)) 


def calculate_baseline(days):
    import S1_methane_equilibrium  as me
    return me.calculate_equilibrium_solubility(days)

def calculate_scenarios(new_depth,days,sc): 

    ''' get sum flux for the scenario
       bubble sizes'''
    scenario_dict = {                  
                     'BS_Basic_seep':                [[4]*4,0.5],  
                     'BS_no_ice':                    [[4]*4,0.5],                    
                     'IOI_Increased_oxidation_rate': [[4]*4,0.5],                    
                     'IMR_Increased_mixing_rate':    [[4]*4,0.5],
                     'SB_Small_bubbles':             [[2]*32,0.5],
                     'RF_Reduced_flux':              [[4]*2,0.5],
                     'IF_Increased_flux':            [[4]*8,0.5]
                     }
    
    df_sum = make_df_sum(s = scenario_dict[sc])

    ## Interpolate to roms depths and 365 days
    f_flow  = interpolate.interp1d(df_sum.depth.values,
                                   df_sum.met_flow.values,
                                   fill_value = 'extrapolate',
                                   kind='nearest') 

    int_flow = f_flow(new_depth)

    smoothed_flow = lowess(new_depth,int_flow,frac=1,is_sorted = False)[:,0]     # here frac is the fraction of the data used when estimating each y-value.
    smoothed_flow_vol,difs_x, volumes  = rate_to_flow(smoothed_flow,new_depth)   # millimole m3 sec 

    # Calculate means,fluxes etc.    
    ma = df_sum.met_cont.max()
    mi = df_sum.met_cont.min()
    to_atm  = np.around(mi,decimals=2)

    base_area = 10# area of the base, m2 
    to_atm_vol = to_atm / (base_area*np.abs(new_depth[1] - new_depth[0])) 
    flux_to_atm = to_atm / base_area #millimole/m2/sec 
    perc = np.around((mi/ma)*100,decimals=2)
    perc1 = np.around(100-perc, decimals=2)
    
    print ('Sum flow to water milliM/sec from the whole! seep {}'.format(sc),df_sum.met_flow.sum())    
    #print ('Mean flow to water milliM/m2/sec from one horizont {}'.format(sc),df_sum.met_flow.mean())
    print ('Scen {} Flux To atmosphere  {} milliM/m2/sec '.format(sc, flux_to_atm), 'max content milliM',ma,'min content in the bubbles milliM',mi)
    print ('CH4 dissolved in the water during uplifting,scenario:{} {} %'.format(sc, perc1)) 
    print ('Percentage of flux CH4 to atm {}  {}%'.format(sc, perc)) 
    print ("     ")    

    table = pd.DataFrame(data = {'depth [m]'          : new_depth,
                                'rate diss [mmol/sec]': smoothed_flow,
                                'difs_x [m]'          : difs_x, 
                                'volumes [m3] '       : volumes ,
                                'flow [mmol/m3 sec]'  : smoothed_flow_vol})

    table = table.sort_values(by = ['depth [m]'])
    table['flow with ice [mmol/m3 sec]'] = table['flow [mmol/m3 sec]']
    table.at[0,'flow with ice [mmol/m3 sec]'] += to_atm_vol  

    word_table(table,sc)
    
    flx_v = table['flow [mmol/m3 sec]'].values
    flx_v_ice = table['flow with ice [mmol/m3 sec]'].values

    flow = np.zeros(shape = (len(new_depth),len(days))) 
    df_flow = pd.DataFrame(flow,index = new_depth,columns = days)
  
    start_wat = 236 # 27 aug
    stop_wat = 303 # 28 oct 

    frac = scenario_dict[sc][1] 
    if frac == 1 :
        for n in days:
            if n < start_wat or n > stop_wat: 
                df_flow[n] = flx_v_ice
            else:
                df_flow[n] = flx_v    
               
    else:     
        import random
        random.seed(14)
        lst = range(1,366)
        n = int(365*frac)
        rand = random.sample(lst,n)
        
        for n in rand:
            if ((n < start_wat) or (n> stop_wat)) and sc != 'BS_no_ice':  
                df_flow[n] = flx_v_ice
            else:
                df_flow[n] = flx_v      

    return df_flow.T #, df_cont.T

if __name__ == '__main__': 


    scenario_dict = {                    
                     'BOM_Basic_seep':  [[4]*4,0.5],  
                     'BOM_no_ice':      [[4]*4,0.5],                  
                     'S_Small_bubbles': [[2]*32,0.5],
                     'F_Reduced_flux':  [[4]*3,0.5],
                     'F2_Reduced_flux': [[4]*2,0.5]}   

    s = scenario_dict['BOM_Basic_seep']
    print (make_df_sum(s).head())
