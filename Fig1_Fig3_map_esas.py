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

from mpl_toolkits.axes_grid1 import make_axes_locatable
# coordinates of needed station 
st_lon = 126.82
st_lat = 76.77

seep_area_lats = [76.5, 77.5, 77.5, 76.5]
seep_area_lons = [121., 121., 132., 132.]
wod_path =  r'Data/Laptev_WOD_area_70_150_f_deeper.nc'

def map_attr(m):
    m.drawmapboundary(fill_color='#cae2f7')
    m.fillcontinents(color='#bdad95')
    m.drawcoastlines(linewidth=0.5)
    m.drawlsmask(ocean_color='#caedfd')

def get_batymnetry():
    etopo1name=r'Data\ETOPO1_Ice_g_gmt4.nc'
    etopo1 = xr.open_dataset(etopo1name)
    #etopo1  = etopo1.where(etopo1.z < 1, drop = True)   
    etopo1  = etopo1.where(((etopo1.lon > 6) & (etopo1.lon < 9)), drop = True)
    etopo1  = etopo1.where(((etopo1.lat > 53.5) & (etopo.lat < 55.2)), drop = True)    
    return etopo1 

def add_bathymetry(etopo1,m,axis,case = 'Laptev'):    


    clevs = np.arange(-300,25,25)                   
    df = etopo1.to_dataframe().unstack()
    lon_e = list(zip(*df.columns))[1]
    lat_e = df.index.values
    z =  df.values 
    
    X_e, Y_e = np.meshgrid(lon_e, lat_e)
    x_e, y_e = m(X_e, Y_e)     
    cs = m.contourf(x_e,y_e,z,clevs,cmap=plt.get_cmap('Blues_r'),
                    extend="both",alpha = 0.9) 

    if case == 'Laptev':   
        m.contour(x_e,y_e,z,clevs,colors = '#225684',
                  linewidths = 0.76,linestyles = 'dotted')
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size="5%", pad = 0.5)
        cb = fig.colorbar(cs,cax=cax,extend='both',ticks=clevs)
        cb.set_label(r"Depth, m", labelpad= 20, y=0.5, rotation=270)
        cax.set_yticklabels(clevs)     

def add_wod(m): 
    ncfile = (wod_path)  
    dss = xr.open_dataset(ncfile)
    lon, lat = m(dss.longitude.values,dss.latitude.values)
    m.scatter(lon,lat,marker='o', alpha = 0.7, s = 5,c = 'k', edgecolors ='k',zorder = 10)
    
def add_polygon(m,axis,a=0.7,c = '#db4832'):
    X, Y = m(seep_area_lons, seep_area_lats)
    poly = Polygon(((X[0],Y[0]),(X[1],Y[1]),(X[2],Y[2]),(X[3],Y[3])),
                   facecolor=c , alpha=a,zorder = 9)
    axis.add_patch(poly) 

def fig1_polar_map(save = False):
    global fig 
    fig = plt.figure(figsize=(4.69,3.27), dpi=100,facecolor='None',edgecolor='None')
    gs = gridspec.GridSpec(1, 1)
    
    gs.update(wspace = 0.07, top = 0.94, bottom = 0.1, 
               hspace = 0,  right = 0.9, left = 0.01) 
    ax_poly = fig.add_subplot(gs[0])

    map_poly = Basemap(projection='npaeqd',boundinglat=70,lon_0=120,resolution='l', ax = ax_poly)
    map_poly.drawparallels(np.arange(70.,99.,5) ) # draw parallels  labels = [1,0,0,0]
    map_poly.drawmeridians(np.arange(0.,365.,20) , labels = [0,1,0,1]) #l r t b 
    map_attr(map_poly)
    etopo = get_batymnetry()
    add_bathymetry(etopo,map_poly,ax_poly,case = 'world')
    add_polygon(map_poly,ax_poly,a=1,c = '#ff3911') 
    
    x_es,y_es = map_poly(155, 75)
    x_l,y_l = map_poly(110, 78)
    x_k,y_k = map_poly(78, 75)
    
    a,f = 0.3,10
    ax_poly.text(x_es, y_es, 'East \nSiberian \nSea',
                fontsize=f,fontweight='bold',
                bbox=dict(facecolor='w', alpha=a ,edgecolor='none'),
                ha='left',va='bottom',color='k')
    ax_poly.text(x_l, y_l, 'Laptev \nSea',fontsize=f,
                fontweight='bold',
                bbox=dict(facecolor='w', alpha=a ,edgecolor='none'), 
                        ha='left',va='bottom',color='k')
    if save == True: 
        plt.savefig('Fig/Figure1_Polar_map.pdf',
                    format =  'pdf',
                    transparent = False)
    else: 
        plt.show()    


def fig1_3_polar_map(save = False):
    global fig  
    etopo = get_batymnetry()

    fig = plt.figure(figsize=(8,3.27), dpi=100,facecolor='None',edgecolor='None')
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2])
    
    gs.update(wspace = 0.03, top = 0.94, bottom = 0.1, 
               hspace = 0,  right = 0.9, left = 0.01) 

    ax_poly = fig.add_subplot(gs[0])
    ax_laptev = fig.add_subplot(gs[1])

    map_poly = Basemap(projection='npaeqd',boundinglat=70,lon_0=120,resolution='l', ax = ax_poly)


    map_laptev = Basemap(resolution='i',projection='poly',\
                lat_0=76,lon_0=125,
                llcrnrlon = 110 ,llcrnrlat = 70,
                urcrnrlon = 150 ,urcrnrlat = 79.2, ax = ax_laptev)

    map_poly.drawparallels(np.arange(70.,99.,5) ) 
    map_poly.drawmeridians(np.arange(0.,180.,20) , labels = [0,0,0,1]) 
    map_attr(map_poly)
    #add_bathymetry(etopo,map_poly,ax_poly,case = 'world')    
    #add_bathymetry(map_poly,ax_poly,case = 'world')
    add_polygon(map_poly,ax_poly,a=0.5,c = '#ff3911') 

    area_lats = [71, 79.2, 79.2, 71]
    area_lons = [109, 109, 150., 150.]

    X, Y = map_poly(area_lons, area_lats)
    poly = Polygon(((X[0],Y[0]),(X[1],Y[1]),(X[2],Y[2]),(X[3],Y[3])),fill = False,
                   edgecolor='k',zorder = 8)
    ax_poly.add_patch(poly)     


    X_lapt, Y_lapt = map_laptev(area_lons, area_lats)    
    
    x_laptev,y_laptev = map_laptev(120,75)
    x1_poly,y1_poly = map_poly(120,75)
    '''for n in np.arange(0,4):
        print (n)
        con = ConnectionPatch(xyA=[X_lapt[n],Y_lapt[n]], xyB=[X[n],Y[n]], coordsA="data", coordsB="data",
                        axesA=ax_laptev, axesB=ax_poly, color="k",zorder = 1)'''
    con = ConnectionPatch(xyA=[300,1.056e6], xyB=[X[1],Y[1]], coordsA="data", coordsB="data",
                    axesA=ax_laptev, axesB=ax_poly, color="k",zorder = 1,alpha = 0.5)       

    con0 = ConnectionPatch(xyA=[X_lapt[0],Y_lapt[0]], xyB=[X[0],Y[0]], coordsA="data", coordsB="data",
                        axesA=ax_laptev, axesB=ax_poly, color="k",zorder = 1) 

    ax_laptev.add_artist(con) 
    ax_laptev.add_artist(con0) 


    map_laptev.drawparallels(np.arange(60.,99.,2) , labels = [0,1,0,0]) 
    map_laptev.drawmeridians(np.arange(0.,365.,10), labels = [0,0,0,1]) 
    map_poly.drawparallels(np.arange(70.,99.,5) ) 
    map_poly.drawmeridians(np.arange(0.,180.,20) , labels = [0,0,0,1]) 

    
    map_attr(map_laptev)
    add_bathymetry(etopo,map_laptev,ax_laptev)
    add_wod(map_laptev)
    add_polygon(map_laptev,ax_laptev) 
    x_l,y_l = map_laptev(120, 74)

    a = 0.3
    f = 10    
    ax_laptev.text(x_l, y_l, 'Laptev Sea',fontsize=f, fontweight='bold',
                bbox=dict(facecolor='w', alpha=a ,edgecolor='none'), 
                        ha='left',va='bottom',color='k')    

    x_es,y_es = map_poly(155, 75)
    x_l,y_l = map_poly(110, 78)
    x_k,y_k = map_poly(78, 75)
    
    a,f = 0.3,10
    ax_poly.text(x_es, y_es, 'East \nSiberian \nSea',
                fontsize=f,fontweight='bold',
                bbox=dict(facecolor='w', alpha=a ,edgecolor='none'),
                ha='left',va='bottom',color='k')
    ax_poly.text(x_l, y_l, 'Laptev \nSea',fontsize=f,
                fontweight='bold',
                bbox=dict(facecolor='w', alpha=a ,edgecolor='none'), 
                        ha='left',va='bottom',color='k')



    if save == True: 
        plt.savefig('Fig/Figure1_Polar_map.pdf',
                    format =  'pdf',
                    transparent = False)
    else: 
        plt.show()         

def fig3_laptev_map(save = False):
    global fig 
    fig = plt.figure(figsize=(4.69,3.27), dpi=100,facecolor='None',edgecolor='None')
    gs = gridspec.GridSpec(1, 1)
    
    gs.update(wspace = 0.07, top = 0.94, bottom = 0.1, 
               hspace = 0,  right = 0.9, left = 0.01) 
    ax_laptev = fig.add_subplot(gs[0])
    
    map_laptev = Basemap(resolution='i',projection='poly',\
                lat_0=76,lon_0=125,
                llcrnrlon = 110 ,llcrnrlat = 70,
                urcrnrlon = 150 ,urcrnrlat = 79.2, ax = ax_laptev)
    map_laptev.drawparallels(np.arange(60.,99.,2) , labels = [0,1,0,0]) 
    map_laptev.drawmeridians(np.arange(0.,365.,10), labels = [0,0,0,1]) 

    map_attr(map_laptev)
    add_bathymetry(map_laptev,ax_laptev)
    add_wod(map_laptev)
    add_polygon(map_laptev,ax_laptev)



    x_l,y_l = map_laptev(120, 74)
    a = 0.3
    f = 10    
    ax_laptev.text(x_l, y_l, 'Laptev Sea',fontsize=f, fontweight='bold',
                bbox=dict(facecolor='w', alpha=a ,edgecolor='none'), 
                        ha='left',va='bottom',color='k')    
  
    if save == True: 
        plt.savefig('Fig/Figure3_Laptev_map.pdf',
                    format =  'pdf',
                    transparent = False)  
    else: 
        plt.show() 


#fig3_laptev_map(True)
#fig1_polar_map(True)
fig1_3_polar_map(False)
#plt.show()