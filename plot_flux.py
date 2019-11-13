import S5_Calculate_scenarios as scen
import pandas as pd 
import numpy as np       
from scipy import interpolate 
from statsmodels.nonparametric.smoothers_lowess import lowess
import xarray as xr




def calculate_scenarios(new_depth,sc): 

    ''' get sum flux for the scenario
       bubble sizes'''
    scenario_dict = {                  
                     'BS_Basic_seep':                [4,4,0.5],  
                     'BS_no_ice':                    [4,4,0.5],                    
                     'IOI_Increased_oxidation_rate': [4,4,0.5],                    
                     'IMR_Increased_mixing_rate':    [4,4,0.5],
                     'SB_Small_bubbles':             [2,32,0.5],
                     'RF_Reduced_flux':              [4,2,0.5],
                     'IF_Increased_flux':            [4,8,0.5],
                     'IF_Increased_flux_2':          [4,4000,0.5],
                     'IA_Increased_activity':        [4,40,1]
                     }
    
    df_sum = scen.make_df_sum(s = scenario_dict[sc])

    ## Interpolate to roms depths and 365 days
    f_flow  = interpolate.interp1d(df_sum.depth.values,
                                   df_sum.met_flow.values,
                                   fill_value = 'extrapolate',
                                   kind='nearest') 

    
    #print ([0,5,10,40,70],new_depth )
    int_flow = f_flow(list(new_depth))
    #print (int_flow)
    smoothed_flow = lowess(new_depth,int_flow,frac=1,is_sorted = False)[:,0]     # here frac is the fraction of the data used when estimating each y-value.
    smoothed_flow_vol,difs_x, volumes  = scen.rate_to_flow(smoothed_flow,new_depth)   # millimole m3 sec 

    # Calculate means,fluxes etc.    
    ma = df_sum.met_cont.max()
    mi = df_sum.met_cont.min()
    to_atm  = np.around(mi,decimals=2)

    base_area = 10 # area of the base, m2 
    to_atm_vol = mi / (base_area*np.abs(new_depth[1] - new_depth[0])) 
    flux_to_SWI = ma / base_area #millimole/m2/sec 
    flux_to_atm = to_atm / base_area #millimole/m2/sec 
    perc = np.around((mi/ma)*100,decimals=2)
    perc1 = np.around(100-perc, decimals=2)
    
    print ('Sum flow to water milliM/sec from the whole! seep {}'.format(sc),df_sum.met_flow.sum())    
    #print ('Mean flow to water milliM/m2/sec from one horizont {}'.format(sc),df_sum.met_flow.mean())
    print ('Scen {} Flux To atmosphere  {} milliM/m2/sec '.format(sc, flux_to_atm), 'max content milliM',ma,'min content in the bubbles milliM',mi)
    print ('Flux to SWI',flux_to_SWI)
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


    
    flx_v = table['flow [mmol/m3 sec]'].values
    flx_v_ice = table['flow with ice [mmol/m3 sec]'].values

    return flx_v,flx_v_ice

ds_baseline = xr.open_dataset(r"C:\Users\ELP\OneDrive\Python_workspace\arctic2030\Data\for-fluxes-plot\Laptev_Baseline_B.nc")     
z =  ds_baseline.depth.values


flx_v,flx_v_ice = calculate_scenarios(z,"BS_Basic_seep")

#basic_scen = xr.open_dataset(r"C:\Users\ELP\OneDrive\Python_workspace\arctic2030\Data\for-fluxes-plot\water_norelax.nc")
basic_scen = xr.open_dataset(r"C:\Users\ELP\OneDrive\Python_workspace\arctic2030\Data\for-fluxes-plot\water.nc")
basic_scen_dic = basic_scen.B_C_DIC[:, 5:]

#print(basic_scen_dic)
dif_dic = basic_scen_dic.values - ds_baseline.DIC.values 
start_wat = 236 # 27 aug
stop_wat = 302 # 28 oct 
import matplotlib.pyplot as plt 
#plt.plot(flx_v,z)

df_gr1 = pd.DataFrame(data = dif_dic[:120,:])     # 1 jan , 30 Apr Winter Mixing 
df_gr2 = pd.DataFrame(data = dif_dic[130:230,:])  # 10 May 18 Aug Stratification 
df_gr3 = pd.DataFrame(data = dif_dic[235:315,:])  # 18 Aug  11 Nov Open water 


import pandas as pd 
step = 3
a = 0.3

plt.fill_betweenx(z,df_gr1.min(axis = 0),df_gr1.max(axis = 0),color = '#469BB9',alpha = a,label = 'Winter mixing')
plt.fill_betweenx(z,df_gr2.min(axis = 0),df_gr2.max(axis = 0),color = '#5AB657',alpha = a,label = "Stratified water")
plt.fill_betweenx(z,df_gr3.min(axis = 0),df_gr3.max(axis = 0),color = '#FF5733',alpha = a,label = "Open water")

plt.plot(df_gr1.mean(axis = 0),z,color = '#469BB9',alpha = 1, linewidth = 5,zorder = 10)
plt.plot(df_gr2.mean(axis = 0),z,color = '#5AB657',alpha = 1, linewidth = 5)
plt.plot(df_gr3.mean(axis = 0),z,color = '#FF5733',alpha = 1, linewidth = 5)

for n in range(0,len(df_gr1.index),step):
    plt.plot(df_gr1.T[n],z,'--',color = '#469BB9',alpha = 0.1,linewidth = 0.5)
for n in range(0,len(df_gr2.index),step):
    plt.plot(df_gr2.T[n],z,'--',color = '#5AB657',alpha = 0.1,linewidth = 0.5)
for n in range(0,len(df_gr3.index),step):
    plt.plot(df_gr3.T[n],z,'--',color = '#FF5733',alpha = 0.1,linewidth = 0.5)

#plt.plot(flx_v_ice,z)
#plt.plot(ds_baseline.DIC.values,z)
plt.title('raise of DIC due to CH4 influx')
plt.ylabel('Depth, m')
plt.xlabel('DIC, $\mu$M')
plt.grid(True)
plt.ylim(80,0)
plt.legend()
plt.show()