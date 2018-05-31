
'''
Created on 25. mai 2018

@author: ELP
'''

import os, time
import datetime
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d,interp2d
import numpy as np
import pandas as pd 
from matplotlib import gridspec


nc_format = 'NETCDF3_CLASSIC' 
f = Dataset('Data\ROMS_Laptev_Sea_{}_east_each_day.nc'.format(nc_format)) 


depth =  f.variables['depth'][:]
depth2 =  f.variables['depth2'][:]

ocean_time =  f.variables['time'][:]
units = f.variables['time'].units
clnd = f.variables['time'].calendar

sal = f.variables['sal'][:][:]
temp = f.variables['temp'][:][:]
kz = f.variables['Kz_s'][:]
hice = f.variables['hice'][:]

fig = plt.figure(figsize=(11,6),dpi = 100)
gs = gridspec.GridSpec(2, 2)
#gs.update(right = 0.99,hspace = 0.5 )
ax1,ax2,ax3,ax4 = plt.subplot(gs[0]),plt.subplot(gs[1]),plt.subplot(gs[2]),plt.subplot(gs[3])
#fig,(ax1,ax2,ax3) = plt.subplots(3)
Times,V_depth = np.meshgrid(ocean_time,depth)
Times2,V_depth2 = np.meshgrid(ocean_time,depth2)


Dates = num2date(Times,units= units,calendar= clnd)
Dates2 = num2date(Times2,units= units,calendar= clnd)
Dates3 = num2date(ocean_time,units= units,calendar= clnd)

start = datetime.datetime(1995, 1, 1, 23, 59)
end = datetime.datetime(1999, 1, 1, 23, 59)

rng = pd.date_range(start, end, freq='A-DEC')
rng2 = pd.date_range(start, end, freq='6M')
#print (start)
#print (rng2)
import matplotlib.dates as mdates

ax1.set_title('Salinity [psu]')
ax2.set_title('Temperature [C]')
ax3.set_title(r'Salinity vertical diffusion coefficient [$meter^2 \cdot second ^{-1}$]')
ax4.set_title(r'ice thickness [m]')
cs1 = ax1.pcolormesh(Dates,V_depth,sal.T,cmap = plt.get_cmap('coolwarm')) #,
cs2 = ax2.pcolormesh(Dates,depth,temp.T,cmap = plt.get_cmap('bwr'),vmin = -5,vmax = 5) 
cs3 = ax3.pcolormesh(Dates2,depth2,kz.T,cmap = plt.get_cmap('Reds')) 


ax4.set_xlim(start,end)
ax4.set_ylim(-0,2.5)

#ax4.axhline(y = 0,c ='w')
ax4.fill_between(Dates3, hice, 3, facecolor='#97C5C3')

#ax4.plot(Dates3,hice,c ='k')
ax4.vlines(rng,0,3,linestyle = '--', lw = 0.5)
ax4.set_xticks(rng)

for axis in (ax1,ax2,ax3):
    axis.set_xticks(rng)
    axis.vlines(rng,80,0,linestyle = '--', lw = 0.5)
    axis.set_xlim(start,end)
    axis.set_ylim(80,0)
   


plt.colorbar(cs1,ax= ax1)    
plt.colorbar(cs2,ax= ax2)   
plt.colorbar(cs3,ax= ax3) 
plt.colorbar(cs3,ax= ax4,ticks = [])

'''
ax1.set_xticks(rng2)
ax1.vlines(rng,80,0,linestyle = '--', lw = 0.5)
ax1.set_xlim(start,end)
ax1.set_ylim(80,0)
           
ax2.set_xticks(rng2)             
ax2.vlines(rng,80,0,linestyle = '--', lw = 0.5)
ax2.set_xlim(start,end)
ax2.set_ylim(80,0)

ax3.set_xticks(rng2)
ax3.vlines(rng,80,0,linestyle = '--', lw = 0.5)
ax3.set_xlim(start,end)
ax3.set_ylim(80,0)'''

#ax2.pcolormesh(ocean_time,depth,temp.T) #,    
#               vmin=0, vmax=0.1) #, alpha = 0.9, label = 'interp')    
plt.tight_layout()
#plt.savefig(r'C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Figures_schemes\tskz_roms_ersem.png')
plt.show()
