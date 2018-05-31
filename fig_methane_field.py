import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline  
from scipy import interpolate 
import methane_equilibrium as me 
import os
from netCDF4 import Dataset
import tkinter as tk # python3
from tkinter.filedialog import askopenfilename # python 3
#To show only the dialog without any other GUI elements
root = tk.Tk()
root.withdraw()

#1  Data from Schroder, Hartwig; Damm, Ellen (2006): 
#Methane concentrations in the water column of the
# Laptev Sea measured on water bottle samples during 
#POLARSTERN cruise ARK-XIV/1b (Transdrift-V). PANGAEA, 

path = r'C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Pangaea_data\ARK-XIV_1b_CH4.tab'
# read table into dataframe 
d = pd.read_table(path, sep='\t',skiprows = 50)
d = d[d['Latitude']<77] 

grouped = d.groupby('Event')
dep_gr = d.groupby('Depth water [m]')
#
stations = grouped.groups.keys()
fig, ax = plt.subplots(figsize = (3,4),dpi = 100)

#lines = []
for st in stations: 
    depth = grouped.get_group(st)['Depth water [m]'] 
    conc =  grouped.get_group(st)['CH4 [nmol/l]'] 
    lines, = ax.plot(conc,depth,'o-', alpha = 0.7,
                     c = '#7f7f7f',  #c4a893',
                     label = 'Polarstern')
    
                                                  
#2 Data to calculate solubility 
file = r'C:\Users\elp\OneDrive\Python_workspace\arctic2030\Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_each_day.nc'        
fh = Dataset(file)     
depth =  fh.variables['depth'][:]

temp =  fh.variables['temp'][40,:] # random day
sal =  fh.variables['sal'][40,:]
fg = 1.8*10**(-6) # 1.87987 ppb to check 
methane2 = []
for n in range(0,len(temp)):
    t = temp[n]  
    s = sal[n]
    d = depth[n]
    met2 = me.calc_methane_depth(t,s,fg,d)[0] 
    methane2.append(met2)

import numpy as np
solubility, = ax.plot(methane2,depth,'ko--',label= r'растворимость')

arr = np.vstack((methane2,depth)).T
np.savetxt('methane_solubility.dat', (arr), delimiter=' ')

#ax.legend(lines[1:], ['Polarstern 1998'],
#          loc='upper right', frameon=False)

#ax.legend(solubility,['solubility']) 
plt.legend(handles = [lines,solubility])   
ax.set_title(r'CH$_4$ [nmol $\cdot$ l$^{-1}$] \n') 
ax.set_ylim(82,0)
plt.show()

#fig_solubility()

'''

for st in stations: 
    d = grouped.get_group(st)['Depth water [m]'] 
    conc =  grouped.get_group(st)['CH4 [nmol/l]'] 
    lines, = ax1.plot(conc,d,'o-', alpha = 0.7,
                     c = 'k',
                     label = 'Polarstern 1998')


solubility, = ax1.plot(methane2,depth,'o--',label= 'Solubility' )
ax1.scatter(d2_ch4,d2_depth,c = '#c4a893')
ax1.set_title('CH$_4$ \n[nmol $\cdot$ l$^{-1}$]') 
ax1.annotate('solubility', xy = (methane2[32],depth[32]))
'''


def fig2():

    
    fig2, (ax2,ax3) = plt.subplots(1, 2, figsize = (9,5))
    path2 = r'C:\Users\elp\OneDrive - NIVA\NIVA_Projects\PERMAFLUX.TRK\Pangaea_data\Bussmann_2016.tab'
    d2 = pd.read_table(path2, sep='\t',skiprows = 46)
    d2 = d2[['CH4 [nmol/l] (in water)','Depth water [m]',
            'k CH4 [1/day]','MOX [nmol/l/day]']]
    
    #print (d2['k CH4 [1/day]'])
    d2 = d2.dropna()
    
    d2 = d2[d2['k CH4 [1/day]'] != '#0.000']
    #print ('after',d2['k CH4 [1/day]'])
    #d2['k CH4 [1/day]'] == 'NaN' 
    d2_ch4 = d2['CH4 [nmol/l] (in water)'].astype('float')
    d2_depth = d2['Depth water [m]'].astype('float')
    d2_k_ch4 = d2['k CH4 [1/day]'].astype('float')
    d2_mox = d2['MOX [nmol/l/day]'].astype('float')
    
    print ('k_ch4', min(d2_k_ch4),max(d2_k_ch4)) 
    for axis in (ax2,ax3):
        axis.set_ylim(20,0)
    
    ax2.set_title('CH$_4$ Fractional turnover \nrate [$day^{-1}$]') 
    ax3.set_title('CH$_4$ oxidation rate \n[nmol $\cdot$ l$^{-1}\cdot day^{-1} $]') 
    
    sc = ax2.scatter(d2_k_ch4,d2_depth,c = d2_ch4)#'#c4a893')
    cax = plt.axes([0.92, 0.1, 0.02, 0.78])
    clb = plt.colorbar(sc,cax = cax)
    clb.set_label('CH$_4$ \n[nmol/l]', labelpad=-30, y=1.1, rotation=0)
    ax2.axvline(d2_k_ch4.median())
    ax2.axvline(d2_k_ch4.mean(),c= 'k')
    ax3.scatter(d2_mox,d2_depth,c = d2_ch4)
    
    
    ax2.annotate('mean', xy = (d2_k_ch4.mean(),10))
    ax2.annotate('median', xy = (d2_k_ch4.median(),15))
    ax3.axvline(0.098,c = 'r',alpha = 0.2)
    ax3.annotate('shakhova', xy = (0.098,10))
    ax3.axvspan(0.098-0.13,0.098+0.13, alpha=0.2, color='red')
    
    ax2.set_xticks([0,0.02,0.04,0.06])
    ax2.set_xlim(0,0.06)
    print ('mean', d2_mox.mean(), 'median', d2_k_ch4.median())
    plt.show()



#plt.show()