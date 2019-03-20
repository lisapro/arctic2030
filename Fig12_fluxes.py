
import matplotlib.dates as mdates
import os
from netCDF4 import Dataset,num2date,date2num,date2index
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
sns.set()

start = 150
stop = 365

fh_b3 = Dataset(r'E:\Users\ELP\Fortran\Ice_model_bubbles\data_spbm_laptev\B3-0.001\water.nc')
fh_b2 = Dataset(r'E:\Users\ELP\Fortran\Ice_model_bubbles\data_spbm_laptev\B2-0.001\water.nc')
fh_b1 = Dataset(r'E:\Users\ELP\Fortran\Ice_model_bubbles\data_spbm_laptev\B1-0.001\water.nc')

fh_r1 = Dataset(r'E:\Users\ELP\Fortran\Ice_model_bubbles\data_spbm_laptev\R1-0.01\water.nc')

fh_b0 = Dataset(r'E:\Users\ELP\Fortran\Ice_model_bubbles\data_spbm_laptev\baseline-0.001\water.nc')

brom_path = (r'E:\Users\ELP\Fortran\Ice_model_bubbles\data_spbm_laptev\baseline-0.01\Laptev_baseline.nc')


base = Dataset(brom_path)
#var = base['B1_50f'][0,-1]*-60*60*24
#var2 = base['B2_50f'][0,-1]*-60*60*24
#var3  = base['B3_50f'][0,-1]*-60*60*24





t = base['time'][start:stop]
t = num2date(t,units = base['time'].units,
                    calendar= 'standard') 

ice_data = Dataset(r'Data\Laptev_average_year_3year.nc')    
ice = ice_data.variables['hice'][start:stop]  

# correct to ice resolution in the model 
new_ice = []
for n in ice:
    if n < 0.05:
        new_ice.append(0)
    else: 
        new_ice.append(n)    
    
ice = new_ice


ice_time = ice_data.variables['time'][start:stop]    
ice_time = num2date(ice_time,units = base['time'].units,
                    calendar= 'standard') 
time = fh_b1.variables['time'][start:stop]    
units = fh_b1.variables['time'].units   
 
#for key,val in enumerate(ice):
#    print (key,ice_time[key],val)
       
format_time = num2date(time,units = units,
                    calendar= 'standard') 

figure = plt.figure(figsize=(7.27, 4.5), dpi=100,
                        facecolor='None',edgecolor='None')  

layer = -1
ch4_b1 = fh_b1.variables['B_CH4_CH4 _flux'][start:stop,layer]
dic_b1 = fh_b1.variables['B_C_DIC   _flux'][start:stop,layer]

ch4_r1 = fh_r1.variables['B_CH4_CH4 _flux'][start:stop,layer]
dic_r1 = fh_r1.variables['B_C_DIC   _flux'][start:stop,layer]

ch4_b2 = fh_b2.variables['B_CH4_CH4 _flux'][start:stop,layer]
dic_b2 = fh_b2.variables['B_C_DIC   _flux'][start:stop,layer]

ch4_b3 = fh_b3.variables['B_CH4_CH4 _flux'][start:stop,layer]
dic_b3 = fh_b3.variables['B_C_DIC   _flux'][start:stop,layer]



gs = gridspec.GridSpec(3, 1,height_ratios=[1, 2,2])
gs.update(left=0.13, right= 0.98,top = 0.96,bottom = 0.1,
                   wspace=0.03, hspace=0.1)  

ax0 = figure.add_subplot(gs[0])  
ax2 = figure.add_subplot(gs[1])
ax3 = figure.add_subplot(gs[2])

ax0.set_xlim(ice_time[0],ice_time[-1])
ax0.stackplot(ice_time,ice,color = '#7B98A8')
ax2.set_xlim(t[0],t[-1])
ax3.set_xlim(t[0],t[-1])    


import numpy as np
ax3.set_ylim(5,-100)
ax2.set_ylim(5,-250)
ax0.set_ylim(0,2.20)

a = 1



def table4_3():
    # Here I calculate the mean flux and range during ice-free period
    tob1 = ch4_b1[ch4_b1 < 0]
    tob3 = ch4_b3[ch4_b3 < 0]
    tor1 = ch4_r1[ch4_r1 < 0]
    tob2 = ch4_b2[ch4_b2 < 0]
    d = 3
    print ('B1',round(sum(tob1)/65000,d),'(',round(min(tob1)/1000,d),round(max(tob1)/1000,d),')',
           '\nB3',round(sum(tob3)/65000,d),'(',round(min(tob3)/1000,d),round(max(tob3)/1000,d),')',
           '\nR1',round(sum(tor1)/65000,d),'(',round(min(tor1)/1000,d),round(max(tor1)/1000,d),')',
           '\nB2',round(sum(tob2)/65000,d),'(',round(min(tob2)/1000,d),round(max(tob2)/1000,d),')')



ax2.plot(format_time,ch4_b1,'r',  label = 'CH4 flux B1',alpha = a)
ax2.plot(format_time,ch4_r1,'g--',label = 'CH4 flux R1',alpha = a)
ax2.plot(format_time,ch4_b2,'m-', label = 'CH4 flux B2',alpha = a) 
ax2.plot(format_time,ch4_b3,'k--',label = 'CH4 flux B3',alpha = a)

x_text,  y_text = 0.05, 0.93
ax2.text(x_text, y_text, 'A)', transform=ax2.transAxes,
                 fontsize=14, fontweight='bold', va='top', ha='right')
ax3.text(x_text, y_text, 'B)', transform=ax3.transAxes,
                 fontsize=14, fontweight='bold', va='top', ha='right')
ax2.legend() 

ax3.xaxis.set_major_formatter(
        mdates.DateFormatter('%b'))   
ax0.xaxis.set_major_formatter(
        mdates.DateFormatter('%b'))   

ax2.set_ylabel('$\mu M\ CH_4\ m ^2 \cdot day ^{-1} \ $')
ax3.set_ylabel('$\mu M\ C\ m ^2 \cdot day ^{-1} \ $')
ax0.set_ylabel('$ice\  $')
ax0.set_yticks([])
ax0.set_xticks([])
ax2.set_xticklabels([])
   

ax3.plot(format_time,dic_b1,'r',  label = 'DIC flux B1',alpha = a)  
ax3.plot(format_time,dic_r1,'g--',  label = 'DIC flux R1',alpha = a) 
ax3.plot(format_time,dic_b2,'m-', label = 'DIC flux B2',alpha = a) 
ax3.plot(format_time,dic_b3,'k--',label = 'DIC flux B3',alpha = a)

ax3.legend() 
#plt.axhline(0)   


fh_r1.close()
fh_b1.close()
fh_b2.close()
fh_b0.close()  
ice_data.close()  

save = True
if save == True:
    plt.savefig('Fig\Figure12.png',format = 'png',transparent = False) 
else:    
    plt.show()   
