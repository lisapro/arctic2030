import pandas as pd
import matplotlib.pyplot as plt
#from scipy.interpolate import UnivariateSpline  
#from scipy import interpolate 
import S1_methane_equilibrium as met 
import os
from netCDF4 import Dataset
import tkinter as tk 
from tkinter.filedialog import askopenfilename #
import numpy as np
root = tk.Tk()
root.withdraw()


def plot_Figure_3(txt = False):
    
    '''
    Function calculates the methane equilibrium 
    solubility and plots it 
    against field data from Schroder, Hartwig; Damm, Ellen (2006): 
    'Methane concentrations in the water column of the
    Laptev Sea measured on water bottle samples during 
    POLARSTERN cruise ARK-XIV/1b (Transdrift-V). PANGAEA, 
    https://doi.pangaea.de/10.1594/PANGAEA.701302 
    '''
    
    path = r'C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Pangaea_data\ARK-XIV_1b_CH4.tab'
    d = pd.read_table(path, sep='\t',skiprows = 50) # read table into dataframe 
    d = d[d['Latitude']<77] 
    
    grouped = d.groupby('Event')
    dep_gr = d.groupby('Depth water [m]')
    stations = grouped.groups.keys()
    fig, ax = plt.subplots(figsize = (3,4),dpi = 100)
    
    for st in stations: 
        depth = grouped.get_group(st)['Depth water [m]'] 
        conc =  grouped.get_group(st)['CH4 [nmol/l]'] 
        lines, = ax.plot(conc,depth,'o-', alpha = 0.7,
                         c = '#7f7f7f',  
                         label = 'Polarstern')
        
                                                 
    #2 Calculate solubility using temerature and salinity from ROMS  
    depth,temp,sal,methane = met.call_met_profile() 
    
    # Plot figure 
    solubility, = ax.plot(methane,depth,'ko--',label= 'Equilibrium \nsolubility')  
    plt.legend(handles = [lines,solubility])   
    ax.set_title(r'CH$_4$ [nmol $\cdot$ l$^{-1}$] \n') 
    ax.set_ylim(82,0)
        
    if txt == True:
        #3 Save txt file with depth and solubility
        arr = np.vstack((methane,depth)).T
        np.savetxt('Data/methane_solubility.dat', (arr), delimiter=' ')   
      
    plt.show()

if __name__ == '__main__':
    
    plot_Figure_3()
    
    
    