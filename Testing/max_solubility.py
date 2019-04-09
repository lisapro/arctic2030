import xarray as xr
import numpy as np
import pandas as pd

def calc_methane_saturation(temp,sal,d):
    Kh_theta = 1.4e-5 * 101.325 # mol/m3 Pa #Convert to M/atm  
    coef = 1900 #[K]
    T = temp + 273.15
    k_sal = 10**(-(sal/58)*0.127)
    press = 1 + d/10
    conc_sat = Kh_theta*1.e6*press*np.exp((coef)*(1/T - 1/298.15))*k_sal #Micromol/l
    return conc_sat

def calculate_saturation_solubility(days):
    roms_path = '..\Data\Laptev_average_year_3year.nc'
    ds = xr.open_dataset(roms_path) 
    sal = ds.sal.values
    temp = ds.temp.values
    depth = ds.depth.values    
    df_satur = pd.DataFrame(index = ds.depth.values,columns = days)  
    for n in days:
        df_satur[n] = [calc_methane_saturation(temp[n][k],sal[n][k],d) for k,d in enumerate(depth)]          
    return df_satur#.T  
days =  [n for n in range(365)]
s = calculate_saturation_solubility(days)
import matplotlib.pyplot as plt
for n in days:
    plt.plot(s[n],s.index)

mmin = int(s[5].min())
mmax = int(s[5].max())
plt.text(x = 10000, y = 20,s = 'min {}\n max {}'.format(mmin,mmax))
plt.xlim(mmin,mmax)    
plt.ylim(80,0)
plt.show()    

