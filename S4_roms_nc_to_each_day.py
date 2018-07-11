'''
Created on 12. sep. 2017

@author: ELP
'''
# script read netcdf file created in read_roms_nc.py script
# write and interpolate the  data to each day 

import os, time
import datetime
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d,interp2d
import numpy as np

f = Dataset(r'C:\Users\ELP\OneDrive\Python_workspace\arctic2030\Data\ROMS_Laptev_Sea_NETCDF3_CLASSIC_east.nc')

# dimensions
ocean_time =  f.variables['time'][:]
depth =  f.variables['depth'][:]
depth2 =  f.variables['depth2'][:]

#read 2d
sal = f.variables['sal'][:][:]
temp = f.variables['temp'][:][:]
kz =  f.variables['Kz_s'][:][:]
rho =  f.variables['rho'][:][:]

DIC = f.variables['DIC'][:][:]
o2 = f.variables['o2'][:][:]
Alk = f.variables['Alk'][:][:]
po4 = f.variables['po4'][:][:]
no3 = f.variables['no3'][:][:]
nh4 = f.variables['no3'][:][:]
Si = f.variables['Si'][:][:]

Z4_c = (f.variables['Z4_c'][:,:]) 
Z5_c = (f.variables['Z5_c'][:,:]) 
Z5_n = (f.variables['Z5_n'][:,:]) 
Z5_p = (f.variables['Z5_p'][:,:]) 
Z6_c = (f.variables['Z6_c'][:,:]) 
Z6_n = (f.variables['Z6_n'][:,:]) 
Z6_p = (f.variables['Z6_p'][:,:]) 

#read 1d
hice =  f.variables['hice'][:]
snow_thick =  f.variables['snow_thick'][:]
tisrf =  f.variables['tisrf'][:]
swrad =  f.variables['swrad'][:]
swradWm2 =  f.variables['swradWm2'][:]
pCO2atm = f.variables['pCO2atm'][:]

#convert time 
times = num2date(ocean_time,units= 'seconds since 1948-01-01 00:00:00' ,
                                   calendar= 'standard')

newtimes =  np.arange(ocean_time[0], ocean_time[-1:], 86400)
newdates = num2date(newtimes,units= 'seconds since 1948-01-01 00:00:00',
                                   calendar= 'standard')

#interpolate 1D variables 
def interpolate_1d_var(var):
    to_interp_var = interp1d(ocean_time, var)
    return  to_interp_var(newtimes)

interp_hice = interpolate_1d_var(hice) 
interp_snow_thick = interpolate_1d_var(snow_thick) 
interp_tisrf = interpolate_1d_var(tisrf) 
interp_swrad = interpolate_1d_var(swrad) 
interp_swradWm2 = interpolate_1d_var(swradWm2) 
interp_pCO2atm = interpolate_1d_var(pCO2atm)

'''
to_interp_hice = interp1d(ocean_time, hice)
to_interp_snow_thick = interp1d(ocean_time, snow_thick)
to_interp_tisrf = interp1d(ocean_time, tisrf)
    
interp_hice = to_interp_hice(newtimes)
interp_snow_thick = to_interp_snow_thick(newtimes)
interp_tisrf = to_interp_tisrf(newtimes)
'''

# interpolate 2D variables 
def interpolate_2d_var(var):
    to_interp_var = interp2d(ocean_time,sorted(depth),var.T,kind = 'linear')
    interp_var = to_interp_var(newtimes,depth)
    #interp_var = interp_var.T 
    return interp_var.T 

interp_sal = interpolate_2d_var(sal)
interp_temp = interpolate_2d_var(temp)
interp_rho = interpolate_2d_var(rho)
interp_DIC = interpolate_2d_var(DIC)
interp_o2 = interpolate_2d_var(o2)
interp_Alk = interpolate_2d_var(Alk)
interp_po4 = interpolate_2d_var(po4)
interp_no3 = interpolate_2d_var(no3)
interp_nh4 = interpolate_2d_var(nh4)
interp_Si = interpolate_2d_var(Si)
interp_Z4_c = interpolate_2d_var(Z4_c)
interp_Z5_c = interpolate_2d_var(Z5_c)
interp_Z5_n = interpolate_2d_var(Z5_n)
interp_Z5_p = interpolate_2d_var(Z5_p)
interp_Z6_c = interpolate_2d_var(Z6_c)
interp_Z6_n = interpolate_2d_var(Z6_n)
interp_Z6_p = interpolate_2d_var(Z6_p)

'''
to_interp_sal = interp2d(ocean_time,depth,sal.T,kind = 'linear')
interp_sal = to_interp_sal(newtimes,depth)
interp_sal = interp_sal.T

to_interp_temp = interp2d(ocean_time,depth,temp.T,kind = 'linear')
interp_temp = to_interp_temp(newtimes,depth)
interp_temp = interp_temp.T

to_interp_rho = interp2d(ocean_time,depth,rho.T,kind = 'linear')
interp_rho = to_interp_rho(newtimes,depth)
interp_rho = interp_rho.T
'''

to_interp_kz = interp2d(ocean_time,depth2,kz.T,kind = 'linear')
interp_kz = to_interp_kz(newtimes,depth2)
interp_kz = interp_kz.T


# plot to print_nc_all_params 1d
def plot_1d():
    plt.plot(newdates,interp_hice,'o',alpha = 0.7, label = 'interp')
    #plt.plot(times,hice,'--',alpha = 0.9, label = 'raw')
    plt.legend()
    plt.show()

def plot_2d(): # to print_nc_all_params 2d
    fig,(ax1,ax2) = plt.subplots(2)
    Times,V_depth = np.meshgrid(ocean_time,depth)
    Times2,V_depth2 = np.meshgrid(newtimes,depth)
    
    Dates2 = num2date(Times2,units= 'seconds since 1948-01-01 00:00:00',
                                       calendar= 'standard') 
    Dates = num2date(Times,units= 'seconds since 1948-01-01 00:00:00',
                                       calendar= 'standard') 
    
    ax1.pcolormesh(newtimes,depth,interp_temp.T) #,
    ax2.pcolormesh(ocean_time,depth,temp.T) #,    
    #               vmin=0, vmax=0.1) #, alpha = 0.9, label = 'interp')    
    #plt.legend()
    plt.show()


#Write data to new netcdf
#### Create nectdf file 
def write_nc():
    # coordinates of needed station 
    st_lon = 126.82
    st_lat = 76.77
    
    dir_name = 'Data'   
    script_dir = os.path.abspath(os.path.dirname(__file__))
    dir_to_save = os.path.abspath(os.path.join(script_dir,dir_name))    
    if not os.path.isdir(dir_to_save):
        os.makedirs(dir_to_save)    
    nc_format = 'NETCDF3_CLASSIC' 
    f1 = Dataset('{}\ROMS_Laptev_Sea_{}_east_each_day.nc'.format(
        dir_to_save,nc_format), mode='w', format= nc_format)    
    f1.description= (
        "lat=%3.2f,lon=%3.2f file from ROMS data interpolated to 1day timedelta"%(
                     st_lat,st_lon))    
    f1.source = 'Elizaveta Protsenko (elp@niva.no)'    
    f1.history = 'Created ' + time.ctime(time.time())    
    f1.createDimension('time', len(newtimes))
    f1.createDimension('z', len(depth))
    f1.createDimension('z2', len(depth2))
    
    v_depth2 = f1.createVariable('depth2','f8',('z2',), zlib= False)
    v_depth2.long_name = "Z-depth matrix for kz, direction up" ;
    v_depth2.units = "meter"
    v_depth2[:] = depth2
    
    v_depth = f1.createVariable('depth','f8',('z',), zlib= False)
    v_depth.long_name = "Z-depth matrix, direction down" ;
    v_depth.units = "meter"
    v_depth[:] = depth
        
    v_time = f1.createVariable('time', 'f8', ('time',), zlib=False)
    v_time.long_name = 'Time in seconds since 1948-01-01 00:00:00'
    v_time.units = 'seconds since 1948-01-01 00:00:00'
    v_time.field = 'time, scalar, series'
    v_time.calendar='standard'
    v_time[:] = newtimes
    
    v_temp=f1.createVariable('temp', 'f8', ('time','z'), zlib=False)
    v_temp.long_name = "Ocean temperature"
    v_temp.units = "degree Celsius"
    v_temp[:,:] = interp_temp
    
    v_sal = f1.createVariable('sal', 'f8', ('time','z'), zlib=False)
    v_sal.long_name = "Time-averaged salinity"
    v_sal.units = "psu"
    v_sal[:,:] = interp_sal
    
    v_kz = f1.createVariable('Kz_s', 'float', ('time','z2'), zlib=False)
    v_kz.long_name = 'Salinity vertical diffusion coefficient'
    v_kz.units = 'meter2 second-1'
    v_kz[:,:] = interp_kz
    
    v_rho = f1.createVariable('rho', 'f8', ('time','z'), zlib=False)
    v_rho.long_name = 'time-averaged density anomaly'
    v_rho.units = 'kilogram meter-1'
    v_rho[:,:] = interp_rho
    
    v_hice = f1.createVariable('hice', 'f8', ('time',), zlib=False)
    v_hice.long_name = 'time-averaged ice thickness in cell'
    v_hice.units = 'meter'
    v_hice[:] = interp_hice
    
    v_snow_thick = f1.createVariable('snow_thick', 'f8', ('time',), zlib=False)
    v_snow_thick.long_name = 'time-averaged thickness of snow cover'
    v_snow_thick.units = 'meter'        
    v_snow_thick[:] = interp_snow_thick
    
    v_tisrf = f1.createVariable('tisrf', 'f8', ('time',), zlib=False)
    v_tisrf.long_name = 'time-averaged temperature of ice surface'
    v_tisrf.units = 'degree Celsius'
    v_tisrf[:] = interp_tisrf
    
    v_DIC = f1.createVariable('DIC', 'f8', ('time','z'), zlib=False)
    v_DIC.long_name = 'time-averaged carbonate/total dissolved inorganic carbon'
    v_DIC.units = 'mmol C/m^3 '
    v_DIC[:] = interp_DIC
    
    v_o2 = f1.createVariable('o2', 'f8', ('time','z'), zlib=False)
    v_o2.long_name = 'time-averaged oxygen/oxygen '
    v_o2.units = 'mmol O_2/m^3'
    v_o2[:] = interp_o2
    
    v_Alk = f1.createVariable('Alk', 'f8', ('time','z'), zlib=False)
    v_Alk.long_name = 'time-averaged carbonate/total alkalinity"'
    v_Alk.units = 'umol/kg'
    v_Alk[:] = interp_Alk    
        
    v_po4 = f1.createVariable('po4', 'f8', ('time','z'), zlib=False)
    v_po4.long_name = 'time-averaged phosphate/phosphorus'
    v_po4.units = 'mmol P/m^3'
    v_po4[:] = interp_po4
        
    v_no3 = f1.createVariable('no3', 'f8', ('time','z'), zlib=False)
    v_no3.long_name = 'time-averaged nitrate/nitrogen'
    v_no3.units = 'mmol N/m^3'
    v_no3[:] = interp_no3
        
    v_nh4 = f1.createVariable('nh4', 'f8', ('time','z'), zlib=False)
    v_nh4.long_name = 'time-averaged ammonium/nitrogen'
    v_nh4.units = 'mmol N/m^3'
    v_nh4[:] = interp_nh4
    
    v_Si = f1.createVariable('Si', 'f8', ('time','z'), zlib=False)
    v_Si.long_name = 'time-averaged silicate/silicate'
    v_Si.units = 'mmol Si/m^3'
    v_Si[:] = interp_Si


    #mesozooplankton
    v_Z4_c = f1.createVariable('Z4_c', 'f8', ('time','z'), zlib=False)
    v_Z4_c.long_name = 'time-averaged mesozooplankton/carbon'
    v_Z4_c.units = 'mmol C/m^3'
    v_Z4_c[:] = interp_Z4_c
    
    
    #microzooplankton
    v_Z5_c = f1.createVariable('Z5_c', 'f8', ('time','z'), zlib=False)
    v_Z5_c.long_name = 'time-averaged microzooplankton/carbon'
    v_Z5_c.units = 'mmol C/m^3'
    v_Z5_c[:] = interp_Z5_c
    
    v_Z5_n = f1.createVariable('Z5_n', 'f8', ('time','z'), zlib=False)
    v_Z5_n.long_name = 'time-averaged microzooplankton/nitrogen'
    v_Z5_n.units = 'mmol N/m^3'
    v_Z5_n[:] = interp_Z5_n
    
    v_Z5_p = f1.createVariable('Z5_p', 'f8', ('time','z'), zlib=False)
    v_Z5_p.long_name = 'time-averaged microzooplankton/phosphorus'
    v_Z5_p.units = 'mmol P/m^3'
    v_Z5_p[:] = interp_Z5_p
    
    #nanoflagellates
    v_Z6_c = f1.createVariable('Z6_c', 'f8', ('time','z'), zlib=False)
    v_Z6_c.long_name = 'time-averaged nanoflagellates/carbon'
    v_Z6_c.units = 'mmol C/m^3'
    v_Z6_c[:] = interp_Z6_c
    
    v_Z6_n = f1.createVariable('Z6_n', 'f8', ('time','z'), zlib=False)
    v_Z6_n.long_name = 'time-averaged nanoflagellates/nitrogen'
    v_Z6_n.units = 'mmol N/m^3'
    v_Z6_n[:] = interp_Z6_n
    
    v_Z6_p = f1.createVariable('Z6_p', 'f8', ('time','z'), zlib=False)
    v_Z6_p.long_name = 'time-averaged nanoflagellates/phosphorus'
    v_Z6_p.units = 'mmol P/m^3'
    v_Z6_p[:] = interp_Z6_p
            
    v_pCO2atm = f1.createVariable('pCO2atm', 'f8', ('time'), zlib=False)
    v_pCO2atm.long_name = 'time-averaged partial pressure of CO2 in air'
    #v_pCO2atm.units = 'uatm'
    v_pCO2atm[:] = interp_pCO2atm
        
    v_swrad = f1.createVariable('swrad', 'f8', ('time'), zlib=False)
    v_swrad.long_name = 'time-averaged solar shortwave radiation flux'
    v_swrad.units = 'watt meter-2'
    v_swrad.negative_value = 'upward flux, cooling'
    v_swrad.positive_value = 'downward flux, heating'
    v_swrad[:] = interp_swrad
     
    v_swradWm2 = f1.createVariable('swradWm2', 'f8', ('time'), zlib=False)
    v_swradWm2.long_name = 'time-averaged solar shortwave radiation flux (under ice)'
    v_swradWm2.units = 'watt meter-2'
    v_swradWm2[:] = interp_swradWm2
    
 
    f1.close()
    
    
#plot_1d()   
plot_2d()
#write_nc()