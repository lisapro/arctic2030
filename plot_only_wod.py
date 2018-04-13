#!/usr/bin/python
# -*- coding: utf-8 -*-
#@author: ELP

from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
    

fh = Dataset(r'E:\Users\ELP\Python plot\arctic2030\Data\WOD-data\data_from_WOD_COLLECTION_(2017-08-21T16-04-41).nc')

fig = plt.figure(figsize=(7,7))

depth =  (fh.variables['var1'][:])
o2 = (fh.variables['var4'][:][:])
o2 = o2*44.6 # convert to mM

no3 = fh.variables['var7'][:][:]
po4 = fh.variables['var5'][:][:]
si = fh.variables['var6'][:][:]
alk = fh.variables['var12'][:][:]
alk = alk * 1000

time =  fh.variables['date_time']
times = num2date(time[:],time.units)
sal = fh.variables['var3'][:][:]
temp = fh.variables['var2'][:][:]
no2 = fh.variables['var8'][:][:]
# pH = fh.variables['var9'][:][:] no pH in roms
# chl = fh.variables['var10'][:][:] no data 
# plankt_biomass = fh.variables['var11'][:][:] no data 

Ntot = fh.variables['var13'][:][:] #no2+no3
pCO2 = fh.variables['var14'][:][:] 

ax = fig.add_subplot(111)
for n in range(0,176):
    ax.plot(o2[n],depth[n],'o-', label = 'wod',c = '#a00a3d') #,s = 10,zorder = 10)  
ax.set_ylim(82,0)
plt.show()