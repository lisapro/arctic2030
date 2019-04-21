import numpy as np
import matplotlib.dates as mdates
import os
from netCDF4 import Dataset,num2date,date2num,date2index
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import xarray as xr
sns.set(style = 'darkgrid')

start = 179
stop = 350

b_path = r'C:\Users\elp\OneDrive - NIVA\BROM_linux_output\with_capping'
brom_path = r'C:\Users\elp\OneDrive - NIVA\BROM_linux_output\Baseline_B\Laptev_Baseline_B.nc'

ice_data =xr.open_dataset(r'Data\Laptev_average_year_3year.nc')
ice = ice_data['hice'][start:stop]
ice = ice.where(ice > 0.05, other = 0)
ice_time = ice_data['time'][start:stop].values

def get_dtst(dname):
    return xr.open_dataset(r'{}\{}\water.nc'.format(b_path,dname))

fh_b =  get_dtst('B-Basic-seep')
fh_s =  get_dtst('S-Small-bubbles')
fh_mi = get_dtst('MI-Increased-horizontal-mixing')
fh_mr = get_dtst('MR-Reduced-horizontal-mixing')
fh_oi = get_dtst('OI-Increased-Oxidation-Rate')
fh_fr = get_dtst('FR-Reduced-flux')

time = fh_b['time'][start:stop].values    
figure = plt.figure(figsize=(7.27, 4.5), dpi=100,
                        facecolor='None',edgecolor='None')  

layer = -1 #surface layer

ch4_name = 'B_CH4_CH4 _flux'
dic_name = 'B_C_DIC   _flux'

def getvar(df):
    return df[ch4_name][start:stop,layer],df[dic_name][start:stop,layer]

ch4_b,dic_b =   getvar(fh_b)
ch4_mi,dic_mi = getvar(fh_mi)
ch4_mr,dic_mr = getvar(fh_mr)
ch4_fr,dic_fr = getvar(fh_fr)
ch4_oi,dic_oi = getvar(fh_oi)
ch4_s,dic_s =   getvar(fh_s)

gs = gridspec.GridSpec(3, 1,height_ratios=[1, 2,2])
gs.update(left=0.13, right= 0.98,top = 0.96,bottom = 0.1,
                   wspace=0.03, hspace=0.1)  

ax0 = figure.add_subplot(gs[0])  
ax1 = figure.add_subplot(gs[1])
ax2 = figure.add_subplot(gs[2])

ax0.set_xlim(ice_time[0],ice_time[-1])
ax0.stackplot(ice_time,ice,color = '#7B98A8')

ax1.set_xlim(time[0],time[-1])
ax2.set_xlim(time[0],time[-1])    

ax2.set_ylim(5,-110)
ax1.set_ylim(5,-410)
ax0.set_ylim(0,2.20)

a = 0.7

x_text,  y_text = 0.05, 0.93
ax1.text(x_text, y_text, 'A)', transform=ax1.transAxes,
                 fontsize=14, fontweight='bold', va='top', ha='right')
ax2.text(x_text, y_text, 'B)', transform=ax2.transAxes,
                 fontsize=14, fontweight='bold', va='top', ha='right')

[axis.xaxis.set_major_formatter(mdates.DateFormatter('%b')) for axis in [ax1,ax2]] 

ax1.set_ylabel('$\mu M\ CH_4\ m ^2 \cdot day ^{-1} \ $')
ax2.set_ylabel('$\mu M\ C\ m ^2 \cdot day ^{-1} \ $')
ax0.set_ylabel('$ice\  $')
ax0.set_yticks([])
ax0.set_xticks([])
ax1.set_xticklabels([])
   
dfs = [dic_b,dic_s,dic_fr,dic_mi,dic_mr,dic_oi]
ch4s = [ch4_b,ch4_s,ch4_fr,ch4_mi,ch4_mr,ch4_oi]
cols = [r'#be1931','r',r'#006218',r'#33a0ce',r'#7b3f8b',r'#e0973e']
nms = ['B','S','FR','MI','MR','OI']   

[ax2.plot(time,d,'-',color = cols[v],label = 'DIC flux {}'.format(nms[v]),
        alpha = a) for v,d in enumerate(dfs)]
[ax1.plot(time,d,'-',color = cols[v],label = 'CH4 flux {}'.format(nms[v]),
        alpha = a) for v,d in enumerate(ch4s)]

ax1.legend() 
ax2.legend() 

ice_data.close()  

save = True
if save == True:
    plt.savefig('Fig\Figure8.png',format = 'png',transparent = False) 
else:    
    plt.show()   
