'''
Created on 26. sep. 2017

@author: ELP
'''
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Polygon


# coordinates of needed station 
st_lon = 126.82
st_lat = 76.77

fig = plt.figure(figsize=(11.69,8.27), dpi=100,facecolor='None',edgecolor='None')
gs = gridspec.GridSpec(1,1) #, width_ratios=[2, 1])

gs.update(left = 0.1,right = 0.91,
          bottom = 0.1, top = 0.95,
          wspace = 0.3 ) 

axins = fig.add_subplot(gs[0])

seep_area_lats = [76.5, 77.5, 77.5, 76.5]
seep_area_lons = [121., 121., 132., 132.]

map2 = Basemap(resolution='i',projection='poly',\
            lat_0=76,lon_0=125,
            llcrnrlon = 110 ,llcrnrlat = 70,
            urcrnrlon = 150 ,urcrnrlat = 78, ax = axins)
map2.drawmapboundary(fill_color='#cae2f7')
map2.fillcontinents(color='#bdad95')
map2.drawcoastlines(linewidth=0.5)
meridians = np.arange(0.,365.,10) #[110,120,130,140,150,160]
parallels = np.arange(60.,99.,2) #[70,75,80,90,95]
map2.drawparallels(parallels, labels = [0,1,0,0]) # draw parallels
map2.drawmeridians(meridians, labels = [0,0,0,1]) 
map2.drawlsmask(ocean_color='#caedfd')

ncfile = (r'C:\Users\elp\OneDrive\Python_workspace\arctic2030\Data\Laptev_WOD_area.nc') 
import xarray as xr
from netCDF4 import Dataset,num2date, date2num
dss = xr.open_dataset(ncfile)

etopo1name=r'C:\Users\elp\OneDrive\Python_workspace\arctic2030\Data\ETOPO1_Ice_g_gmt4.nc'
etopo1 = xr.open_dataset(etopo1name)
etopo1  = etopo1.where(etopo1.lon > 100, drop = True)
etopo1  = etopo1.where(etopo1.lon < 150, drop = True)
etopo1  = etopo1.where(etopo1.lat > 70, drop = True)
etopo1  = etopo1.where(etopo1.z < 1, drop = True)
#etopo1  = etopo1.where(etopo1.z > -2000, drop = True)
df = etopo1.to_dataframe()

df = df.unstack()
lon_e = list(zip(*df.columns))[1]
print (len(lon_e))
lat_e = df.index.values
z =  df.values 
print (len(lat_e))
print (z.shape)

X_e, Y_e = np.meshgrid(lon_e, lat_e)
x_e, y_e = map2(X_e, Y_e) 
import matplotlib.colors as colors
clevs = np.arange(-200,-10,5)
clevs2 =[-200,-50]
import cmocean.cm as cmo
#print (np.amin(z),np.amax(z))
#bounds = np.linspace(-500, 0, 100)
#norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

cs = map2.contourf(x_e,y_e,z,clevs,cmap=plt.get_cmap('Blues_r'), extend="both",alpha = 0.8) 
map2.contour(x_e,y_e,z,clevs2,colors = '#225684', linewidths = 0.76,linestyles = 'solid') #plt.get_cmap('pink_r')
fig.colorbar(cs,ax=axins, extend='both')

X, Y = map2( seep_area_lons, seep_area_lats )
lon, lat = map2(dss.longitude.values,dss.latitude.values)
poly = Polygon(((X[0],Y[0]),(X[1],Y[1]),(X[2],Y[2]),(X[3],Y[3])),
               facecolor='#db4832', alpha=0.5,zorder = 1 )
axins.add_patch(poly)
#axins.text(X[0],Y[0]-200000, r'Laptev sea', fontsize=14)
plt.gca().add_patch(poly)
plt.xticks(visible=False)
plt.yticks(visible=False)

map2.scatter(lon,lat,marker='o', alpha = 0.7, s = 4, color='k',zorder = 10)

plt.savefig('seep_map.png', transparent = True)
#plt.show()