'''
Created on 25. sep. 2017

@author: ELP
'''

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes,inset_axes

from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Polygon


fig = plt.figure(figsize=(6.69,4.27), dpi=100,facecolor='None',edgecolor='None')
gs = gridspec.GridSpec(1,1) #, width_ratios=[2, 1])

gs.update(left = 0.1,right = 0.91,
          bottom = 0.1, top = 0.95,
          wspace = 0.3 ) 

ax = fig.add_subplot(gs[0])
#ax2 = fig.add_subplot(gs[1])

map = Basemap(
            resolution='l',projection='stere',\
            lat_0=76,lon_0=125,
            llcrnrlon = 120 ,llcrnrlat = 75.9,
            urcrnrlon = 140 ,urcrnrlat = 78,ax = ax)

#map = Basemap(projection='npstere',boundinglat=70,lon_0=120,resolution='l',ax = ax)

map.drawmapboundary(fill_color='#cae2f7')
map.fillcontinents(color='#bdad95')
map.drawcoastlines(linewidth=0.5)
map.drawlsmask(ocean_color='#caedfd')
meridians = np.arange(0.,365.,10) #[110,120,130,140,150,160]
parallels = np.arange(70.,99.,2) #[70,75,80,90,95]
map.drawparallels(parallels, labels = [1,0,0,0]) # draw parallels
map.drawmeridians(meridians, labels = [0,0,0,1]) 
map.drawmapscale(134,78.1, 134, 78.1,  200, barstyle='fancy')
'''
map2.drawmapboundary(fill_color='#cae2f7')
map2.fillcontinents(color='#bdad95')
map2.drawcoastlines(linewidth=0.5)
lat=76.74,lon=126.71 file from ROMS data for station (lat=76.77,lon=126.82)
'''
weight = 'semibold'

# coordinates of needed station 
st_lon = 126.82
st_lat = 76.77

roms_lon = 126.71
roms_lat = 76.74

roms2_lon = 127.0
roms2_lat = 76.56


seep_x,seep_y = map(st_lon,st_lat)
roms2_x,roms2_y = map(roms2_lon,roms2_lat)
#roms2_x,roms2_y = map(roms2_lon,2)

cax = map.scatter(st_lon, st_lat,latlon= True , color ='#880303',
              label = "Porewater data (Nutrients) \nfrom seep station", 
          edgecolor = '#49361f', s= 75, zorder = 5, alpha = 0.8)

cax2 = map.scatter(roms_lon, roms_lat,latlon= True , color ='y',
              label = "Input data from ROMS \n(T,S,Kz,Ice,Snow)", 
          edgecolor = '#49361f', s= 75, zorder = 10, alpha = 0.8)

cax2 = map.scatter(roms2_lon, roms2_lat,latlon= True , color ='y',
              label = "Input data from ROMS \n(T,S,Kz,Ice,Snow)", 
          edgecolor = '#49361f', s= 75, zorder = 10, alpha = 0.8)

lena_delta_x, lena_delta_y = map(125,73)

# Add patch 
# patch = patches.PathPatch(path, facecolor='#b7ebd9', lw=1,alpha = 0.1)
# ax.add_patch(patch)
seep_area_lats = [76.5, 77.,77.5, 77.5, 76.5]
seep_area_lons = [121., 121.,130., 132., 132.]
x, y = map( seep_area_lons, seep_area_lats )

x = np.array(x)
y = np.array(y)
#xy = np.concatenate((x,y), axis=0)
#xy = np.ravel(xy)
#xy = xy.reshape(4,2)
#xy = np.array(zip(x,y)).T
poly = Polygon(((x[0],y[0]),(x[1],y[1]),(x[2],y[2]),(x[3],y[3]),(x[4],y[4])), facecolor='#fff7e6', alpha=1,zorder = 1,label = 'Methane seepage rate' )
ax.add_patch(poly)
plt.gca().add_patch(poly)
plt.legend(handles=[cax,cax2,poly],loc= 2)

ax.annotate('Seep polygon \n(Shakhova \net.al 2017)', (x[1],y[1]), xytext=(0, -50),weight = weight, fontsize=14, #size= txtsize,
            #arrowprops=dict(arrowstyle="-",connectionstyle="arc3"),
            textcoords='offset points')

ax.annotate('Seep st. \n76.77N \n126.82E ', (seep_x,seep_y), xytext=(35,0) ,weight = weight, fontsize=14, #size= txtsize,
            arrowprops=dict(arrowstyle="-",connectionstyle="arc3"),
            textcoords='offset points')

ax.annotate('ROMS st.\n76.74N \n126.71E ', (roms2_x,roms2_y), xytext=(15, -49),weight = weight, fontsize=14, #size= txtsize,
            arrowprops=dict(arrowstyle="-",connectionstyle="arc3"),
            textcoords='offset points')

'''
axins = inset_axes(ax,width="50%",height="50%") #,bbox_to_anchor=[120,79,5,5] )
axins.set_xlim(121,132)
axins.set_ylim(76.5,77.5)

map2 = Basemap(resolution='i',projection='stere',\
            lat_0=76,lon_0=125,
            llcrnrlon = 110 ,llcrnrlat = 72,
            urcrnrlon = 150 ,urcrnrlat = 78, ax = axins)
map2.drawmapboundary(fill_color='#cae2f7')
map2.fillcontinents(color='#bdad95')
map2.drawcoastlines(linewidth=0.5)
meridians = np.arange(0.,365.,10) #[110,120,130,140,150,160]
parallels = np.arange(70.,99.,2) #[70,75,80,90,95]
map2.drawparallels(parallels, labels = [0,1,0,0]) # draw parallels
map2.drawmeridians(meridians, labels = [0,0,0,1]) 
map2.drawlsmask(ocean_color='#caedfd')
map2.scatter(st_lon, st_lat,latlon= True , color ='#880303',
              label = "seep station", 
          edgecolor = '#49361f', s= 75, zorder = 5, alpha = 0.8)

map2.scatter(roms_lon, roms_lat,latlon= True , color ='y',
              label = "roms station", 
          edgecolor = '#49361f', s= 75, zorder = 10, alpha = 0.8)

X, Y = map2( seep_area_lons, seep_area_lats )
poly = Polygon(((X[0],Y[0]),(X[1],Y[1]),(X[2],Y[2]),(X[3],Y[3])),
                facecolor='#fff7e6', alpha=1,zorder = 1 )
axins.add_patch(poly)
axins.text(X[0],Y[0]-200000, r'Laptev sea', fontsize=14)
plt.gca().add_patch(poly)'''
plt.xticks(visible=False)
plt.yticks(visible=False)
#ax.annotate('Lena delta', (lena_delta_x, lena_delta_y), xytext=(0, 0),weight = weight, #size= txtsize,
#            #arrowprops=dict(arrowstyle="-",connectionstyle="arc3"),
#            textcoords='offset points')
#map.scatter(x,y)

#plt.savefig('laptev_map_corrected_polygon.png', format='png',transparent=True)
#mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
plt.show()