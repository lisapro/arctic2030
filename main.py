'''
Created on 4. sep. 2017

@author: ELP
'''
# Script to read all the nc files in a folder 
# with different time periods from ROMS
# extracts data from one station 
# combines all time steps and writes to new netcdf 
# z recalculated to meters from sigmas 

import os, time
from netCDF4 import MFDataset,Dataset #,num2date,date2num
import numpy as np
import matplotlib.pyplot as plt 
import warnings 

test = False

if test == True: 
    files = ["X:/trondkr/CloudStation/ARCTIC2030/a20_avg_11705_arctic2030.nc",
             "X:/trondkr/CloudStation/ARCTIC2030/a20_avg_11733_arctic2030.nc",
             "X:/trondkr/CloudStation/ARCTIC2030/a20_avg_11761_arctic2030.nc"]
    f = MFDataset(files)
else:
    f = MFDataset("X:/trondkr/CloudStation/ARCTIC2030/*.nc")



latitude = np.array(f.variables['lat_rho'])
longitude = np.array(f.variables['lon_rho'])

# coordinates of needed station 
st_lon = 126.82
st_lat = 76.77


# function 'def find_xi_eta' is based on 
# Model2roms  Python toolbox     
# https://github.com/trondkr/model2roms 

## finds xi and eta of the station closest to the needed position
def find_xi_eta(st_lon,st_lat):
    distance = np.zeros((longitude.shape),dtype=np.float64)
    listd=[]
    
    '''First, create a list of distances from the station of interest, while '''
    ''' also save the matrix of distances that contains the info to get the '''
    ''' index pair that the distance of interest corresponds to'''
    
    for eta in range(len(latitude[:,0])):
        for xi in range(len(latitude[0,:])):
            distance[eta,xi] = np.sqrt(
                (latitude[eta,xi]-st_lat)**2.0 + (longitude[eta, xi] - st_lon)**2.0 )
            listd.append(distance[eta,xi])     
    value = min(listd)
    itemindex = np.where(distance==value)
    return itemindex

itemindex = find_xi_eta(st_lon,st_lat)
eta = itemindex[0] 
xi = itemindex[1]

# Read the data from needed station 
temp = np.array(f.variables['temp'][:,:,eta,xi])
temp = temp[:,:,0,0]
sal = np.array(f.variables['salt'][:,:,eta,xi])
sal = sal[:,:,0,0]
ntime = np.array(f.variables['ocean_time'][:])

sigma = np.array(f.variables['s_rho'][:])

kz = np.array(f.variables['AKs'][:,:,eta,xi]) 
kz = kz[:,:,0,0]

rho = np.array(f.variables['rho'][:,:,eta,xi]) 
rho = rho[:,:,0,0]

hice = np.array(f.variables['hice'][:,eta,xi]) 
hice = hice[:,0,0]
snow_thick = np.array(f.variables['snow_thick'][:,eta,xi]) 
snow_thick = snow_thick[:,0,0]
tisrf = np.array(f.variables['tisrf'][:,eta,xi]) 
tisrf = tisrf[:,0,0]

#long_name = "salinity vertical diffusion coefficient"
       

# Read the variables needed to recalculate 
# depth from sigma values  to meters

vtransform = f.variables['Vtransform'][:]
vstretching = f.variables['Vstretching'][:]
h = f.variables["h"][eta,xi]
hc = f.variables["hc"][:]
Tcline = hc # 250.0  
theta_s = f.variables['theta_s'][:]
theta_b = f.variables['theta_b'][:]
Nlevels = 40

f.close()
  



# this part of the script is taken from 
# Model2roms  Python toolbox     
# https://github.com/trondkr/model2roms

class s_coordinate(object):
    """
    Song and Haidvogel (1994) vertical coordinate transformation (Vtransform=1) and
    stretching functions (Vstretching=1).
    
    return an object that can be indexed to return depths

    s = s_coordinate(h, theta_b, theta_s, Tcline, N)
    """

    def __init__(self, h, theta_b, theta_s, Tcline, N, vtransform, vstretching, zeta=None):
      
        self.h = np.asarray(h)
        self.hmin = h.min()
        self.theta_b = theta_b
        self.theta_s = theta_s
        self.Tcline = Tcline
        self.N = int(N)
        self.Np = self.N+1
        self.vtransform = vtransform
        self.vstretching = vstretching
        self.hc = min(self.hmin, self.Tcline)

        self.Vtrans = 1

        if self.vtransform==1:
            if (self.Tcline > self.hmin):
                warnings.warn(
                    'Vertical transformation parameters are not defined correctly ',
                    'in either gridid.txt or in the history files: \n Tcline = %d and ', 
                    'hmin = %d. \n You need to make sure that Tcline <= hmin when using ',
                    'transformation 1.' %(self.Tcline,self.hmin))

        self.c1 = 1.0
        self.c2 = 2.0
        self.p5 = 0.5       

        if zeta is None:
            self.zeta = np.zeros(h.shape)
        else:
            self.zeta = zeta
        
        self._get_s_rho()
        self._get_s_w()
        self._get_Cs_r()
        self._get_Cs_w()

        self.z_r = z_r(self.h, self.hc, self.N, self.s_rho, self.Cs_r, self.zeta, self.Vtrans)
        self.z_w = z_w(self.h, self.hc, self.Np, self.s_w, self.Cs_w, self.zeta, self.Vtrans)


    def _get_s_rho(self):
        lev = np.arange(1,self.N+1,1)
        ds = 1.0 / self.N
        self.s_rho = -self.c1 + (lev - self.p5) * ds

    def _get_s_w(self):
        lev = np.arange(0,self.Np,1)
        ds = 1.0 / (self.Np-1)
        self.s_w = -self.c1 + lev * ds

    def _get_Cs_r(self):
        if (self.theta_s >= 0):
            Ptheta = np.sinh(self.theta_s * self.s_rho) / np.sinh(self.theta_s)
            Rtheta = np.tanh(self.theta_s * (self.s_rho + self.p5)) / \
                      (self.c2 * np.tanh(self.p5 * self.theta_s)) - self.p5
            self.Cs_r = (self.c1 - self.theta_b) * Ptheta + self.theta_b * Rtheta
        else:
            self.Cs_r = self.s_rho

    def _get_Cs_w(self):
        if (self.theta_s >= 0):
            Ptheta = np.sinh(self.theta_s * self.s_w) / np.sinh(self.theta_s)
            Rtheta = np.tanh(self.theta_s * (self.s_w + self.p5)) / \
                      (self.c2 * np.tanh(self.p5 * self.theta_s)) - self.p5
            self.Cs_w = (self.c1 - self.theta_b) * Ptheta + self.theta_b * Rtheta
        else:
            self.Cs_w = self.s_w
            
class s_coordinate_4(s_coordinate):
    """
    A. Shchepetkin (2005) UCLA-ROMS vertical coordinate transformation (Vtransform=2) and
    stretching functions (Vstretching=4).
    
    return an object that can be indexed to return depths

    s = s_coordinate_4(h, theta_b, theta_s, Tcline, N)
    """

    def __init__(self, h, theta_b, theta_s, Tcline, N,  vtransform, vstretching, zeta=None):
        self.h = np.asarray(h)
        self.hmin = h.min()
        self.theta_b = theta_b
        self.theta_s = theta_s
        self.Tcline = Tcline
        self.N = int(N)
        self.Np = self.N+1
        self.vtransform = vtransform
        self.vstretching = vstretching

        self.hc = self.Tcline

        self.Vtrans = 4

        self.c1 = 1.0
        self.c2 = 2.0
        self.p5 = 0.5

        if zeta is None:
            self.zeta = np.zeros(h.shape)
        else:
            self.zeta = zeta

        self._get_s_rho()
        self._get_s_w()
        self._get_Cs_r()
        self._get_Cs_w()

        self.z_r = z_r(self.h, self.hc, self.N, self.s_rho, self.Cs_r, self.zeta, self.Vtrans)
        self.z_w = z_w(self.h, self.hc, self.Np, self.s_w, self.Cs_w, self.zeta, self.Vtrans)
        

    def _get_s_rho(self):
        super(s_coordinate_4, self)._get_s_rho()

    def _get_s_w(self):
        super(s_coordinate_4, self)._get_s_w()

    def _get_Cs_r(self):
        if (self.theta_s > 0):
            Csur = (self.c1 - np.cosh(self.theta_s * self.s_rho)) / \
                     (np.cosh(self.theta_s) - self.c1)
        else:
            Csur = -self.s_rho**2
        if (self.theta_b > 0):
            Cbot = (np.exp(self.theta_b * Csur) - self.c1 ) / \
                   (self.c1 - np.exp(-self.theta_b))
            self.Cs_r = Cbot
        else:
            self.Cs_r = Csur         

    def _get_Cs_w(self):
        if (self.theta_s > 0):
            Csur = (self.c1 - np.cosh(self.theta_s * self.s_w)) / \
                     (np.cosh(self.theta_s) - self.c1)
        else:
            Csur = -self.s_w**2
        if (self.theta_b > 0):
            Cbot = (np.exp(self.theta_b * Csur) - self.c1 ) / \
                   ( self.c1 - np.exp(-self.theta_b) )
            self.Cs_w = Cbot
        else:
            self.Cs_w = Csur

class z_r(object):
    """
    return an object that can be indexed to return depths of rho point

    z_r = z_r(h, hc, N, s_rho, Cs_r, zeta, Vtrans)
    """

    def __init__(self, h, hc, N, s_rho, Cs_r, zeta, Vtrans):
        self.h = h
        self.hc = hc
        self.N = N
        self.s_rho = s_rho
        self.Cs_r = Cs_r
        self.zeta = zeta
        self.Vtrans = Vtrans

    def __getitem__(self, key):

        if isinstance(key, tuple) and len(self.zeta.shape) > len(self.h.shape):
            zeta = self.zeta[key[0]]
            res_index = (slice(None),) + key[1:]
        elif len(self.zeta.shape) > len(self.h.shape):
            zeta = self.zeta[key]
            res_index = slice(None)
        else:
            zeta = self.zeta
            res_index = key
     
        if self.h.ndim == zeta.ndim:       # Assure a time-dimension exists
            zeta = zeta[np.newaxis, :]
        
        ti = zeta.shape[0]
        z_r = np.empty((ti, self.N) + self.h.shape, 'd')
        if self.Vtrans == 1:
            for n in range(ti):
                for  k in range(self.N):
                    z0 = self.hc * self.s_rho[k] + (self.h - self.hc) * self.Cs_r[k]
                    z_r[n,k,:] = z0 + zeta[n,:] * (1.0 + z0 / self.h)
        elif self.Vtrans == 2 or self.Vtrans == 4:
            for n in range(ti):
                for  k in range(self.N):
                    z0 = (self.hc * self.s_rho[k] + self.h * self.Cs_r[k]) / \
                          (self.hc + self.h)
                    z_r[n,k,:] = zeta[n,:] + (zeta[n,:] + self.h) * z0

        return np.squeeze(z_r[res_index])
    
class z_w(object):
    """
    return an object that can be indexed to return depths of w point

    z_w = z_w(h, hc, Np, s_w, Cs_w, zeta, Vtrans)
    """

    def __init__(self, h, hc, Np, s_w, Cs_w, zeta, Vtrans):
        self.h = h
        self.hc = hc
        self.Np = Np
        self.s_w = s_w
        self.Cs_w = Cs_w
        self.zeta = zeta
        self.Vtrans = Vtrans

    def __getitem__(self, key):

        if isinstance(key, tuple) and len(self.zeta.shape) > len(self.h.shape):
            zeta = self.zeta[key[0]]
            res_index = (slice(None),) + key[1:]
        elif len(self.zeta.shape) > len(self.h.shape):
            zeta = self.zeta[key]
            res_index = slice(None)
        else:
            zeta = self.zeta
            res_index = key
     
        if self.h.ndim == zeta.ndim:       # Assure a time-dimension exists
            zeta = zeta[np.newaxis, :]
      
        ti = zeta.shape[0]
        z_w = np.empty((ti, self.Np) + self.h.shape, 'd')
        if self.Vtrans == 1:
            for n in range(ti):
                for  k in range(self.Np):
                    z0 = self.hc * self.s_w[k] + (self.h - self.hc) * self.Cs_w[k]
                    z_w[n,k,:] = z0 + zeta[n,:] * (1.0 + z0 / self.h)
        elif self.Vtrans == 2 or self.Vtrans == 4:
            for n in range(ti):
                for  k in range(self.Np):
                    z0 = (self.hc * self.s_w[k] + self.h * self.Cs_w[k]) / \
                          (self.hc + self.h)
                    z_w[n,k,:] = zeta[n,:] + (zeta[n,:] + self.h) * z0

        return np.squeeze(z_w[res_index])    

vgrid = s_coordinate_4(h, theta_b, theta_s, Tcline, Nlevels, vtransform, vstretching, zeta=None)
depth = - vgrid.z_r[0,:] # depth in m (axis from surface to bottom)
depth2 = - vgrid.z_w[0,:] # interfaces depths
#z_w = z_w
#z_r = vgrid.z_r[0,:]
#Cs_rho = vgrid.Cs_r
#self.s_w = vgrid.s_w
#Cs_w = vgrid.Cs_w
#s_rho = vgrid.s_rho





### plot to check  
#for n in range(0,12):
#    plt.plot(temp[n],z_r)
#plt.ylim(200,0)
#plt.show()


#### Create nectdf file 
nc_format = 'NETCDF3_CLASSIC' 
f1 = Dataset('ROMS_Laptev_Sea_{}.nc'.format(nc_format), mode='w', format= nc_format)

f1.description=" lat=%3.2f,lon=%3.2f file from ROMS data"%(st_lat,st_lon)
f1.source = 'Elizaveta Protsenko (elp@niva.no)'

f1.history = 'Created ' + time.ctime(time.time())

f1.createDimension('time', len(ntime))
f1.createDimension('z', len(depth))
f1.createDimension('z2', len(depth2))

v_depth = f1.createVariable('depth','f8',('z',), zlib= False)
v_depth.long_name = "Z-depth matrix, direction down" ;
v_depth.units = "meter"
v_depth[:] = depth

v_depth2 = f1.createVariable('depth2','f8',('z2',), zlib= False)
v_depth2.long_name = "Z-depth matrix for kz, direction down" ;
v_depth2.units = "meter"
v_depth2[:] = depth2

v_time = f1.createVariable('time', 'f8', ('time',), zlib=False)
v_time.long_name = 'seconds since 1948-01-01 00:00:00'
v_time.units = 'seconds since 1948-01-01 00:00:00'
v_time.field = 'time, scalar, series'
v_time.calendar='standard'
v_time[:] = ntime

v_temp=f1.createVariable('temp', 'f8', ('time','z'), zlib=False)
v_temp.long_name = "Ocean temperature"
v_temp.units = "degree Celsius"
v_temp[:,:] = temp

v_sal = f1.createVariable('sal', 'f8', ('time','z'), zlib=False)
v_sal.long_name = "Time-averaged salinity"
v_sal.units = "psu"
v_sal[:,:] = sal

v_kz = f1.createVariable('Kz_s', 'f8', ('time','z2'), zlib=False)
v_kz.long_name = 'Salinity vertical diffusion coefficient'
v_kz.units = 'meter2 second-1'
v_kz[:,:] = kz

v_rho = f1.createVariable('rho', 'f8', ('time','z'), zlib=False)
v_rho.long_name = 'time-averaged density anomaly'
v_rho.units = 'kilogram meter-1'
v_rho[:,:] = rho

v_hice = f1.createVariable('hice', 'f8', ('time',), zlib=False)
v_hice.long_name = 'time-averaged ice thickness in cell'
v_hice.units = 'meter'
v_hice[:] = hice

v_snow_thick = f1.createVariable('snow_thick', 'f8', ('time',), zlib=False)
v_snow_thick.long_name = 'time-averaged thickness of snow cover'
v_snow_thick.units = 'meter'
v_snow_thick[:] = hice

v_tisrf = f1.createVariable('tisrf', 'f8', ('time',), zlib=False)
v_tisrf.long_name = 'time-averaged temperature of ice surface'
v_tisrf.units = 'degree Celsius'
v_tisrf[:] = tisrf

f1.close()


if __name__ == '__main__':
    pass