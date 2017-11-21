# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 14:45:45 2017

@author: zagoven
"""
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import matplotlib.pyplot as plt
import numpy as np 
import datetime 
import matplotlib.gridspec as gridspec

plt.style.use('ggplot')
#choose the year,month and day in our netcdf
year=2014
month=11
day=24
dt=datetime.datetime(year,month,day,12,0)
# input netcdf:
filein='ROMS_Laptev_Sea_NETCDF3_CLASSIC_each_day.nc'

# output file for using in methan.f90 model:
fileout='Laptev_sea_{0}_{1}_{2}'.format(day,month, year)+'.dat'

# output temperature plot:
plot_T='Temp_Laptev_sea_{0}_{1}_{2}'.format(day,month, year)+'.png' 
# output salinity plot:
plot_S='Sal_Laptev_sea_{0}_{1}_{2}'.format(day,month, year)+'.png'
# create the human formate for our time array
dt=datetime.datetime(year,month,day,12,0)

#reading the netcdf
f = Dataset(filein)
depth =  f.variables['depth'][:] # means all value
sal = f.variables['sal'][:]
temp = f.variables['temp'][:]
ocean_time =  f.variables['time'][:]
times = num2date(ocean_time,
                 units= 'seconds since 1948-01-01 00:00:00' ,
                                   calendar= 'standard')

#np.where(times==dt) #ищет в массивах, где наше время совпадает с dt
#sort the depth for methan.f90 
#(not sure that nessesary - need to check)
True_depth = -np.sort(depth)
ind=np.where(times==dt)[0][0] #index in time, [0] for the acsess to the 3 element

sal[ind] #salinity from one particular time moment
temp[ind] # the same for temperature

def save_datfile():
    #create 3 columns from 3 arrays
    out_array=np.column_stack((True_depth,sal[ind],temp[ind])) 
    np.savetxt(fileout,out_array, delimiter=' ') #save the main data file

# function to plot to subplots in on figure 
def plot_2var(var1,var2):
    figure = plt.figure(figsize=(8,6), dpi=100)
    gs = gridspec.GridSpec(1, 2)
    gs.update(wspace=0.2,hspace = 0.3,left=0.1,
       right=0.97,bottom = 0.05, top = 0.85) 
    
    ax00 = figure.add_subplot(gs[0])
    ax01 = figure.add_subplot(gs[1])    
    ax00.plot(var1,depth)
    ax01.plot(var2,depth)
    
    figure.suptitle('Date: {}'.format(times[ind]), fontsize=14,weight = 'roman')
    ax00.set_title('Temperature', y=1.06, fontsize=13)
    #y is a position of title 1 is upper border of picture 
    ax01.set_title('Salinity', y=1.06, fontsize=13)
    
    for axis in (ax00,ax01) :
        axis.set_ylim(150,0)
        axis.xaxis.tick_top() 
        axis.set_ylabel('Depth,m') 
    plt.show()
    
    
def plot_var(var,title,filename):
    plt.plot(var,depth) #creation temperature plot
    plt.ylim(150,0)
    plt.title(str(times[ind])+title)
    plt.xlabel(title +'°С')
    plt.ylabel('Depth,m')
    #it is nessesary to put .savefig before .show
    plt.savefig(filename, dpi=400, facecolor='w', edgecolor='r',
            orientation='portrait', papertype=None, forma=None,
            transparent=False, bbox_inches=None, pad_inches=0.1,
            frameon=None)
    plt.show()


#Uncomment to call functions 

#plot_var(temp[ind],' Temperature',plot_T) 
#plot_var(sal[ind],' Salinity',plot_S)
#save_datfile()

plot_2var(temp[ind],sal[ind])

