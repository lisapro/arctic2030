'''
Created on 28. jun. 2017

Module for calculation the differences between 
chosen scenario and baseline scenario 
Used for creation the table 4 in the paper

@author: ELP
'''

import os,sys,datetime 
from netCDF4 import Dataset,num2date,date2num,date2index
import numpy as np
import numpy.ma as ma
from dateutil.relativedelta import relativedelta

def read_var(name):
    var_water = np.array(
        fh_water.variables[name][:]).T  
    var_water_base  = np.array(
        fh_water_base.variables[name][:]).T     
    var_water_dif = var_water-var_water_base              
    data_units = fh_water.variables[name].units                  
    return var_water_dif #,data_units

def make_dif(water_fname,water_base_fname,case):      

    global fh_water,fh_water_base
    fh_water_base = Dataset(water_base_fname)

    fh_water =  Dataset(water_fname)   

    depth_water = np.array(fh_water.variables['z'][:])    
    min_water = np.amin(depth_water)
    max_water = np.amax(depth_water)
            

    time = fh_water.variables['time']      
    time2 = fh_water.variables['time'][:]
    time_units = fh_water.variables['time'].units
    format_time = num2date(time2,units = time_units,calendar= 'standard')

    start_year = 1991
    stop_year = 1992 

    to_start = datetime.datetime(start_year,1,1,12,0)
    to_stop= datetime.datetime(stop_year,1,1,12,0)
            
    start = date2index(to_start,time,
                        calendar=None, select='nearest')
    stop = date2index(to_stop, time,
                        calendar=None, select='nearest') 

    def to_count(vname):
        var_water_dif = read_var(vname)
        v = var_water_dif[:,start:stop]
        t = time2[start:stop]
        from numpy import unravel_index

        def get_td(ind):
            d = depth_water[ind[0]]
            tt = t[ind[1]] 
            tt = num2date(tt,units = time_units,calendar= 'standard') 
            return d,tt

        indices_max = unravel_index(v.argmax(), v.shape)       
        indices_min = unravel_index(v.argmin(), v.shape)

        d_max,t_max = get_td(indices_max)
        d_min,t_min = get_td(indices_min)

        mm = np.max(v)
        mmin = np.min(v)

        file.write('\nvar: {}, \nmaxdif: {}  \ndepth {}, time {}'.format(vname,mm, d_max,t_max))
        file.write('\nvar: {}, \nmindif: {}  \ndepth {}, time {}'.format(vname,mmin,d_min,t_min))

    file.write('\ncase: {}'.format(case)) 
    [to_count(var) for var in ('B_CH4_CH4','B_BIO_O2_rel_sat','B_C_DIC','B_pH_pH')]    

    fh_water.close()
    
base = r'C:\Users\elp\OneDrive - NIVA\BROM_linux_output'
#base_nocap = r'{}\no_capping'.format(base)
base_cap = r'{}\with_capping'.format(base)

water_base_b_fname = r'{}\baseline-B\water.nc'.format(base)
#water_base_b_fname = r'{}\{}\water.nc'.format(base_cap,'baseline-B-new') 
water_base_o_fname = r'{}\Baseline_O\water.nc'.format(base)
water_B  = r'{}\{}\water.nc'.format(base,'B-Basic-seep')      
water_FR = r'{}\{}\water.nc'.format(base_cap,'FR-Reduced-flux')
water_FR2 = r'{}\{}\water.nc'.format(base_cap,'FR-Reduced-flux')
water_MI = r'{}\{}\water.nc'.format(base_cap,'MI-Increased-horizontal-mixing')       
water_MR = r'{}\{}\water.nc'.format(base_cap,'MR-Reduced-horizontal-mixing')
water_OI = r'{}\{}\water.nc'.format(base_cap,'OI-Increased-Oxidation-Rate')
water_S  = r'{}\{}\water.nc'.format(base_cap,'S-Small-bubbles')      

if __name__ == '__main__':
    file = open("differences.txt", "w") 
    make_dif(water_B, water_base_b_fname,case = 'B')        
    '''make_dif(water_FR2,water_base_b_fname,case = 'FR2')   ,
    make_dif(water_FR,water_base_b_fname,case = 'FR')           
    make_dif(water_MR,water_base_b_fname,case = 'MR')  
    make_dif(water_MI,water_base_b_fname,case = 'MI')     
    make_dif(water_S, water_base_b_fname,case = 'S')          
    make_dif(water_OI,water_base_o_fname,case = 'OI')'''
    file.close() 