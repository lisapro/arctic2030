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
            llcrnrlon = 90 ,llcrnrlat = 65,
            urcrnrlon = 180 ,urcrnrlat = 78, ax = axins)
'''
m = Basemap(llcrnrlon=-35.,llcrnrlat=-30,urcrnrlon=80.,urcrnrlat=50.,\
            resolution='l',area_thresh=1000.,projection='poly',\
            lat_0=0.,lon_0=20.)
'''
map2.drawmapboundary(fill_color='#cae2f7')
map2.fillcontinents(color='#bdad95')
map2.drawcoastlines(linewidth=0.5)
meridians = np.arange(0.,365.,10) #[110,120,130,140,150,160]
parallels = np.arange(60.,99.,2) #[70,75,80,90,95]
map2.drawparallels(parallels, labels = [0,1,0,0]) # draw parallels
map2.drawmeridians(meridians, labels = [0,0,0,1]) 
map2.drawlsmask(ocean_color='#caedfd')



X, Y = map2( seep_area_lons, seep_area_lats )
poly = Polygon(((X[0],Y[0]),(X[1],Y[1]),(X[2],Y[2]),(X[3],Y[3])),
               facecolor='#fff7e6', alpha=1,zorder = 1 )
axins.add_patch(poly)
#axins.text(X[0],Y[0]-200000, r'Laptev sea', fontsize=14)
plt.gca().add_patch(poly)
plt.xticks(visible=False)
plt.yticks(visible=False)

#plt.savefig('small_seep_map.eps', transparent = True)
plt.show()