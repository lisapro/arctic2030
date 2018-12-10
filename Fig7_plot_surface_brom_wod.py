#!/usr/bin/python
# -*- coding: utf-8 -*-
#@author: ELP

from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy.ma as ma

sns.set()
ncfile = r'Data/Laptev_WOD_area_70_150_f_deeper.nc'

fig  = plt.figure(figsize = (7,6))

fh = Dataset(ncfile)
fig.suptitle('Surface concentrations from BROM vs WOD')
time =  fh.variables['date_time']
times = num2date(time[:],time.units)
months = [ time.timetuple()[1] for time in times ]  

depth =  (fh.variables['var1'][:])
o2 = (fh.variables['var4'][:,0])*44.6 # convert to mM
no3 = fh.variables['var7'][:,0]
po4 = fh.variables['var5'][:,0]
si = fh.variables['var6'][:,0]
alk = fh.variables['var12'][:,0]* 1000
sal = fh.variables['var3'][:,0]
temp = fh.variables['var2'][:,0]
pH = fh.variables['var9'][:,0]


alk = ma.masked_less(alk, 1750)
sal = ma.masked_less(sal, 25)
si = ma.masked_greater(si, 50)
po4 = ma.masked_greater(po4, 2)
temp  = ma.masked_greater(temp, 3)

titles = ['$O_2\ \mu M$','$NO_3\ \mu M$','$PO_4\ \mu M$','$Si\ \mu M$',
          '$Alkalinity \ \mu M$','$Temperature \ ^\circ C$','$Salinity\ psu$','$pH$']

vars = [o2,no3,po4,si,alk,temp,sal,pH]
vars_br = ['B_BIO_O2','B_NUT_NO3','B_NUT_PO4','B_NUT_Si',
           'B_C_Alk','temp','sal','B_pH_pH']

dss = xr.open_dataset('Data\water.nc')
group = dss.groupby(dss.time.dt.month)
dss  = group.mean()
dss_min  = group.min()
dss_max  = group.max()
m = dss.month

axx = []
cols = 2
gs = gridspec.GridSpec(len(titles) // cols + 1, cols)
gs.update(hspace=0.9,bottom = -0.05)

for n,title in enumerate(titles):
    row = (n // cols)
    col = n % cols    
    axx.append(fig.add_subplot(gs[row, col]))
    axx[-1].set_title(titles[n])    
    axx[-1].scatter(months,vars[n],alpha = 0.3,c = '#7f7f7f')
    axx[-1].fill_between(m, dss_min[vars_br[n]], dss_max[vars_br[n]],alpha = 0.3,facecolor ='#de7b5c')
    axx[-1].plot(m,dss[vars_br[n]],'o--',alpha = 0.7,c ='#de7b5c') 
    
    axx[-1].set_xticks(np.arange(1,13))
    axx[-1].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    
#plt.show()
plt.savefig('Fig/Figure7_SurfaceBROM_WOD.pdf',
                    format =  'pdf',
                    transparent = False)