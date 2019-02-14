
import matplotlib.dates as mdates
import os
from netCDF4 import Dataset,num2date,date2num,date2index
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
sns.set()
stop = 366


fh_b2 = Dataset(r'E:\Users\ELP\Fortran\Ice_model_bubbles\data_spbm_laptev\B2_50\water.nc')
fh_b1 = Dataset(r'E:\Users\ELP\Fortran\Ice_model_bubbles\data_spbm_laptev\B1_50\water.nc')
fh_b0 = Dataset(r'E:\Users\ELP\Fortran\Ice_model_bubbles\data_spbm_laptev\baseline\water.nc')
brom_path = (r'E:\Users\ELP\Fortran\Ice_model_bubbles\data_spbm_laptev\B1_50\Laptev_baseline.nc')
ice_data = Dataset(r'Data\Laptev_average_year_3year.nc')    
ice = ice_data.variables['hice'][:stop]  
new_ice = []   
for val in ice:
    if val < 0.05:
        val = 0  
        new_ice.append(val)
    else:
        new_ice.append(val)    
ice = new_ice        
for key,val in enumerate(ice):
    print (key,val*100)


ice_time = ice_data.variables['time'][:stop]    

time = fh_b1.variables['time'][:stop]    
units = fh_b1.variables['time'].units   
 

format_time = num2date(time,units = units,
                    calendar= 'standard') 
#form_ice_time = num2date(ice_time,units = units)
figure = plt.figure(figsize=(7.27, 7.1), dpi=100,
                        facecolor='None',edgecolor='None')  

layer = 0
ch4_b1 = fh_b1.variables['B_CH4_CH4 _flux'][:stop,layer]
ch4_b2 = fh_b2.variables['B_CH4_CH4 _flux'][:stop,layer]
ch4_b0_1 = fh_b0.variables['B_CH4_CH4 _flux'][:stop,layer]
ch4_b0 = fh_b0.variables['B_CH4_CH4 _flux'][:stop,0]
ch4_b0_2 = fh_b0.variables['B_CH4_CH4 _flux'][:stop,2]

gs = gridspec.GridSpec(4, 1)
gs.update(left=0.1, right= 0.9,top = 0.96,bottom = 0.1,
                   wspace=0.03, hspace=0.35)  

ax0 = figure.add_subplot(gs[0])
ax1 = figure.add_subplot(gs[1])   
ax2 = figure.add_subplot(gs[2])
ax3 = figure.add_subplot(gs[3])

for ax in (ax0,ax1,ax2):
    ax.set_xlim(format_time[0],format_time[-1])
ax0.stackplot(format_time,ice,color = '#7B98A8')

#ax1.plot(format_time,ch4_b0)
ax1.plot(format_time,ch4_b0_1)
#ax1.plot(format_time,ch4_b0_2)

import numpy as np

#ax2.set_ylim(-10,10)

ax2.plot(format_time,ch4_b1)   
ax2.plot(format_time,ch4_b2) 
ax2.xaxis.set_major_formatter(
        mdates.DateFormatter('%b'))   
   
plt.axhline(0)   

#scen = fh_scen.variables['B2_50f'][:stop,0]

import xarray as xr
ds = xr.open_dataset(brom_path)
f = ds['B2_50f'].loc['1991-1':'1991-12'] 
f = f.where(f != 0)
cm = plt.get_cmap('Greys')
cs = f.plot(x='time', y = 'depth', ax = ax3,
                add_colorbar = False,levels = 50, 
            cmap = cm,vmin = 0,vmax = 0.0003) 
    
    
#t = fh_scen.variables['days'][:stop]
   
#ax3.set_ylim(0,2.30689e-05)
ax3.set_xlim()


fh_b1.close()
fh_b2.close()
fh_b0.close()  
ice_data.close()  
 





 
 
 
 
 
plt.show()   
