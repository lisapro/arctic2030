from datetime import time  
import os, time
import xarray as xr
from netCDF4 import Dataset,num2date, date2num
#import S3_roms_nc_seasonal_mean as make_nc
import S5_Calculate_scenarios as scen 
import numpy as np
import pandas as pd
import datetime

'''
Script to make nc file 
for relaxation brom model scenarios 
based on baseline run of brom model 

BROM reads ROMS and adds BBL to it,
'''

global time_units
time_units = 'seconds since 1948-01-01 00:00:00'
 
start = datetime.datetime(1991,1,1)
sec_start = date2num(start, 
                units = time_units,
                calendar = 'standard')

def d2n(delta): 
    d  = datetime.timedelta(days=int(delta))  
    s = date2num((d+start), 
                units = time_units,
                calendar = 'standard')
    return s 

def make_nc_baseline():
    f_brom = Dataset(r'C:\Users\elp\OneDrive - NIVA\BROM_linux_output\water.nc', mode='r')
    f_roms = Dataset('Data\Laptev_average_year_3year.nc', mode='r')
    
    # remove bbl! 
    z = f_brom.variables['z'][5:]
    z2 = f_brom.variables['z_faces'][5:]
    
    depth = f_roms.dimensions['depth']
    depth2 = f_roms.dimensions['depth2']
    times = f_roms.dimensions['time']
    
    n_years = 3
    
    seconds = [d2n(delta) for delta in np.arange(1,365*n_years+1)]    
    days = [day for day in np.arange(1,365*n_years+1)]

    nc_format = 'NETCDF3_CLASSIC'
    f1 = Dataset('Data\Laptev_baseline.nc', mode='w', format= nc_format)
    f1.description="Baseline after BROM run" 
    f1.source = 'Elizaveta Protsenko (elp@niva.no)'
    f1.history = 'Created ' + time.ctime(time.time())
    
    
    f1.createDimension('time',  size=len(times))
    f1.createDimension('days',  size=len(days))    
    f1.createDimension('depth', size= len(depth))
    f1.createDimension('depth2', size=len(depth2))
    
    v_depth = f1.createVariable('depth','f8',('depth',), zlib= False)
    v_depth.long_name = "Z-depth matrix, direction up" 
    v_depth.units = "meter"
    v_depth[:] = z
    
    v_depth2 = f1.createVariable('depth2','f8',('depth2',), zlib= False)
    v_depth2.long_name = "Z-depth2 matrix, direction up" 
    v_depth2.units = "meter"
    v_depth2[:] = z2

    v_time = f1.createVariable('time', 'f8', ('time',), zlib=False)
    v_time.long_name = 'Time in seconds since 1948-01-01 00:00:00'
    v_time.units = time_units
    v_time.field = 'time, scalar, series'
    v_time.calendar='standard'
    v_time[:] = seconds
    
    v_days = f1.createVariable('days', 'f8', ('days',), zlib=False)
    v_days.long_name = 'number of day in a year'
    v_days.field = 'time, scalar, series'
    v_days[:] = days

    v_o2 = f1.createVariable('O2', 'f8', ('time','depth'), zlib=False)
    v_o2.long_name = 'time-averaged oxygen/oxygen '
    v_o2.units = 'mmol O_2/m^3'
    v_o2[:] = f_brom.variables['B_BIO_O2'][:,5:]   
    
    v_temp = f1.createVariable('temp', 'f8', ('time','depth'), zlib=False)
    v_temp.long_name = 'time-averaged ocean temperatue '
    v_temp.units = 'Celsius'
    v_temp[:] = f_brom.variables['temp'][:,5:]

    v_sal = f1.createVariable('sal', 'f8', ('time','depth'), zlib=False)
    v_sal.long_name = 'time-averaged salinity '
    v_sal.units = 'psu'
    v_sal[:] = f_brom.variables['sal'][:,5:]

    v_rho = f1.createVariable('rho', 'f8', ('time','depth'), zlib=False)
    v_rho.long_name = 'time-averaged density anomaly '
    v_rho.units = 'kilogram meter-3'
    v_rho[:] = f_brom.variables['rho'][:,5:]

    v_Kz_s = f1.createVariable('Kz_s', 'f8', ('time','depth2'), zlib=False)
    v_Kz_s.long_name = 'Salinity vertical diffusion coefficient'
    v_Kz_s.units = 'm^2/sec'
    v_Kz_s[:] = f_brom.variables['Kz_s'][:,5:]
    
    v_hice = f1.createVariable('hice', 'f8', ('time',), zlib=False)
    v_hice.long_name = 'time-averaged ice thickness in cell'
    v_hice.units = 'meter'
    v_hice[:] = f_roms.variables['hice'][:] 

    v_tisrf = f1.createVariable('tisrf', 'f8', ('time',), zlib=False)
    v_tisrf.long_name = 'time-averaged temperature of ice surface'
    v_tisrf.units = 'Celsius'
    v_tisrf[:] = f_roms.variables['tisrf'][:]

    v_swrad = f1.createVariable('swrad', 'f8', ('time',), zlib=False)
    v_swrad.long_name = 'time-averaged solar shortwave radiation flux'
    v_swrad.units = 'watt meter-2'
    v_swrad.negative_value = 'upward flux,cooling'
    v_swrad[:] = f_roms.variables['swrad'][:]

    v_po4 = f1.createVariable('PO4', 'f8', ('time','depth'), zlib=False)
    v_po4.long_name = ' phosphate/phosphorus'
    v_po4.units = 'mmol P/m^3'
    v_po4[:] = f_brom.variables['B_NUT_PO4'][:,5:]
       
    v_no3 = f1.createVariable('NO3', 'f8', ('time','depth'), zlib=False)
    v_no3.long_name = ' nitrate/nitrogen'
    v_no3.units = 'mmol N/m^3'
    v_no3[:] = f_brom.variables['B_NUT_NO3'][:,5:]

    v_no2 = f1.createVariable('NO2', 'f8', ('time','depth'), zlib=False)
    v_no2.long_name = ' nitrite/nitrogen'
    v_no2.units = 'mmol N/m^3'
    v_no2[:] = f_brom.variables['B_NUT_NO2'][:,5:]
    
    v_Si = f1.createVariable('Si', 'f8', ('time','depth'), zlib=False)
    v_Si.long_name = ' silicate/silicate'
    v_Si.units = 'mmol Si/m^3'
    v_Si[:] = f_brom.variables['B_NUT_Si'][:,5:] 
    
    v_dic = f1.createVariable('DIC', 'f8', ('time','depth'), zlib=False)
    v_dic.long_name = 'Dissolved inorganic Carbon '
    v_dic.units = 'mmol c/m^3'
    v_dic[:] = f_brom.variables['B_C_DIC'][:,5:]     

    v_alk = f1.createVariable('Alk', 'f8', ('time','depth'), zlib=False)
    v_alk.long_name = 'Dissolved inorganic Carbon '
    v_alk.units = 'mmol c/m^3'
    v_alk[:] = f_brom.variables['B_C_Alk'][:,5:]  

    v_ch4 = f1.createVariable('CH4', 'f8', ('time','depth'), zlib=False)
    v_ch4.long_name = 'Methane'
    v_ch4.units = 'mmol c/m^3'
    v_ch4[:] = f_brom.variables['B_CH4_CH4'][:,5:] 
    
    import S5_Calculate_scenarios as scen 
    stop = 365*n_years
    days_1 = np.arange(1,367)
    #zeros = scen.calculate_spin_up(z,days_1) 

    def three_years(flux):    
        flux = (pd.concat([flux,flux,flux],
               axis = 0,ignore_index=True)).iloc[:stop,:]       
        return flux     

    #flux_B1_50,cont_B1_50 = scen.calculate_scenarios(z,False,days_1,'B1_50')
    #flux_B2_50,cont_B2_50 = scen.calculate_scenarios(z,False,days_1,'B2_50')
    #flux_B3_50,cont_B3_50 = scen.calculate_scenarios(z,False,days_1,'B3_50')   

    # 1 step create scenario for one year  
    flux_BOM = scen.calculate_scenarios(z,days_1,'BOM_Basic_seep')      
    flux_S = scen.calculate_scenarios(z,days_1,'S_Small_bubbles')                                      
    flux_F = scen.calculate_scenarios(z,days_1,'F_Reduced_flux')                        
                     
    # 2 step write it in three years 
    flux_BOM = three_years(flux_BOM)

    flux_S = three_years(flux_S)
    flux_F = three_years(flux_F)

    v_BOM = f1.createVariable('Scenario_BOM_Basic_seep_flux', 'f8', ('time','depth'), zlib=False)
    v_BOM.long_name = 'Methane inflow scenario B_Basic_seep,O_Increased_oxidation_rate,M_Increased_mixing_rate, 4 bubbles 4mm 50% of time'
    v_BOM.units = 'mmol CH4/m^2 sec'
    v_BOM[:] = flux_BOM

    v_S = f1.createVariable('Scenario_S_flux', 'f8', ('time','depth'), zlib=False)
    v_S.long_name = 'Methane inflow scenario S_Small_bubbles 32 bubbles 2 mm 50% of time'
    v_S.units = 'mmol CH4/m^2 sec'
    v_S[:] = flux_S 

    v_F = f1.createVariable('Scenario_F_flux', 'f8', ('time','depth'), zlib=False)
    v_F.long_name = 'Methane inflow scenario F_Reduced_flux 2 bubbles 4 mm 50% of time'
    v_F.units = 'mmol CH4/m^2 sec'
    v_F[:] = flux_F

    import S1_methane_equilibrium  as me
    slb_year = me.calculate_equilibrium_solubility(days_1) 
    slb_year = three_years(slb_year)

    v_B1_slb = f1.createVariable('CH4_Slb', 'f8', ('time','depth'), zlib=False)
    v_B1_slb.long_name = 'Methane equilibrium solubility'
    v_B1_slb.units = 'mmol/l CH4 '
    v_B1_slb[:]  = slb_year  

    sat_year = me.calculate_saturation_solubility(days_1) 
    sat_year = three_years(sat_year) 
    v_sat = f1.createVariable('Sat_sol', 'f8', ('time','depth'), zlib=False)
    v_sat.long_name = 'Methane saturation solubility'
    v_sat.units = 'micrommol/l CH4 '
    v_sat[:]  = sat_year 


    f1.close()
    f_brom.close()
    f_roms.close()
    

if __name__=='__main__':
    make_nc_baseline()