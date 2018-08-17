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
import xarray as xr
from netCDF4 import Dataset,num2date, date2num
import matplotlib.colors as colors
import cmocean.cm as cmo

# coordinates of needed station 
st_lon = 126.82
st_lat = 76.77

seep_area_lats = [76.5, 77.5, 77.5, 76.5]
seep_area_lons = [121., 121., 132., 132.]

def map_attr(m):
    m.drawmapboundary(fill_color='#cae2f7')
    m.fillcontinents(color='#bdad95')
    m.drawcoastlines(linewidth=0.5)
    m.drawlsmask(ocean_color='#caedfd')

def add_bathymetry(m,axis,case = 'Laptev'):    
    etopo1name=r'C:\Users\elp\OneDrive\Python_workspace\arctic2030\Data\ETOPO1_Ice_g_gmt4.nc'
    etopo1 = xr.open_dataset(etopo1name)
    etopo1  = etopo1.where(etopo1.z < 1, drop = True)   
     
    if case == 'Laptev':
        etopo1  = etopo1.where(etopo1.lon > 100, drop = True)
        etopo1  = etopo1.where(etopo1.lon < 150, drop = True)
        etopo1  = etopo1.where(etopo1.lat > 70, drop = True)
        clevs = np.arange(-300,25,25) 
    else: 
        etopo1  = etopo1.where(etopo1.lat > 62, drop = True)
        clevs = np.arange(-300,100,100) 

                      
    df = etopo1.to_dataframe().unstack()
    lon_e = list(zip(*df.columns))[1]
    lat_e = df.index.values
    z =  df.values 
    
    X_e, Y_e = np.meshgrid(lon_e, lat_e)
    x_e, y_e = m(X_e, Y_e)     
    cs = m.contourf(x_e,y_e,z,clevs,cmap=plt.get_cmap('Blues_r'),
                    extend="both",alpha = 0.9) 
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(axis)
    cax = divider.append_axes("right", size="5%", pad = 0.7)
    cb = fig.colorbar(cs,cax=cax,extend='both',ticks=clevs)
    cb.set_label(r"Depth, m", labelpad= 20, y=0.5, rotation=270)
    cax.set_yticklabels(clevs)     
    if case == 'Laptev':   
        m.contour(x_e,y_e,z,clevs,colors = '#225684',
                  linewidths = 0.76,linestyles = 'dotted')

def add_wod(m): 
    ncfile = (r'C:\Users\elp\OneDrive\Python_workspace\arctic2030\Data\Laptev_WOD_area_70_150_f_deeper.nc')  
    dss = xr.open_dataset(ncfile)
    lon, lat = m(dss.longitude.values,dss.latitude.values)
    m.scatter(lon,lat,marker='o', alpha = 0.7, s = 5,c = 'k', edgecolors ='k',zorder = 10)
    
def add_polygon(m,axis,a=0.7,c = '#db4832'):
    X, Y = m(seep_area_lons, seep_area_lats)
    poly = Polygon(((X[0],Y[0]),(X[1],Y[1]),(X[2],Y[2]),(X[3],Y[3])),
                   facecolor=c , alpha=a,zorder = 9)
    axis.add_patch(poly) 

def fig1_polar_map():
    global fig 
    fig = plt.figure(figsize=(4.69,3.27), dpi=100,facecolor='None',edgecolor='None')
    gs = gridspec.GridSpec(1, 1)
    
    gs.update(wspace = 0.07, top = 0.94, bottom = 0.1, 
               hspace = 0,  right = 0.9, left = 0.01) 
    axins2 = fig.add_subplot(gs[0])
    map = Basemap(projection='npaeqd',boundinglat=70,lon_0=120,resolution='l', ax = axins2)
    map.drawparallels(np.arange(70.,99.,5) ) # draw parallels  labels = [1,0,0,0]
    map.drawmeridians(np.arange(0.,365.,20) , labels = [0,1,0,1]) #l r t b 
    map_attr(map)
    add_bathymetry(map,axins2,case = 'world')
    add_polygon(map,axins2,a=1,c = '#ff3911') 
    
    x_es,y_es = map(155, 75)
    x_l,y_l = map(110, 78)
    x_k,y_k = map(78, 75)
    
    a,f = 0.3,10
    axins2.text(x_es, y_es, 'East \nSiberian \nSea',
                fontsize=f,fontweight='bold',
                bbox=dict(facecolor='w', alpha=a ,edgecolor='none'),
                ha='left',va='bottom',color='k')
    axins2.text(x_l, y_l, 'Laptev \nSea',fontsize=f,
                fontweight='bold',
                bbox=dict(facecolor='w', alpha=a ,edgecolor='none'), 
                        ha='left',va='bottom',color='k')
    plt.savefig('Maps/Figure1_Polar_map.png', transparent = True)
    
def fig3_laptev_map():
    global fig 
    fig = plt.figure(figsize=(4.69,3.27), dpi=100,facecolor='None',edgecolor='None')
    gs = gridspec.GridSpec(1, 1)
    
    gs.update(wspace = 0.07, top = 0.94, bottom = 0.1, 
               hspace = 0,  right = 0.9, left = 0.01) 
    axins = fig.add_subplot(gs[0])
    map2 = Basemap(resolution='i',projection='poly',\
                lat_0=76,lon_0=125,
                llcrnrlon = 110 ,llcrnrlat = 70,
                urcrnrlon = 150 ,urcrnrlat = 79.2, ax = axins)
    map2.drawparallels(np.arange(60.,99.,2) , labels = [0,1,0,0]) 
    map2.drawmeridians(np.arange(0.,365.,10), labels = [0,0,0,1]) 

    map_attr(map2)
    add_bathymetry(map2,axins)
    add_wod(map2)
    add_polygon(map2,axins)   
    x_l,y_l = map2(120, 74)
    a = 0.3
    f = 10    
    axins.text(x_l, y_l, 'Laptev Sea',fontsize=f, fontweight='bold',
                bbox=dict(facecolor='w', alpha=a ,edgecolor='none'), 
                        ha='left',va='bottom',color='k')    
    plt.savefig('Maps/Fig3_Figure_Laptev_map.png', transparent = True)



fig3_laptev_map()
#fig1_polar_map()
#plt.show()