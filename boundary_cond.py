'''
Created on 28. nov. 2017

@author: ELP
'''

from datetime import date
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import numpy.ma as ma
from numpy.ma.core import masked_equal
import xarray as xr

file = r'C:\Users\ELP\workspace\arctic2030\src\Data\data_from_WOD_COLLECTION_(2017-08-21T16-04-41).nc'
data = xr.open_dataset(file)
o2 = data['var4']
print (o2)

def use_netcdf4():
    interpolate = False #True
    clean_masked = False #True
    
    fh = Dataset(file)
    o2 = fh.variables['var4'][:,0]*44.6
    no3 = fh.variables['var7'][:,0]
    po4 = fh.variables['var5'][:,0]
    si = fh.variables['var6'][:,0]
    
    depth = fh.variables['var1'][:,0]
    time = fh.variables['date_time']
    time =  num2date(time[:],time.units)
    #month = time.date.month
    
    day_in_year = []
    for n,item in enumerate(time):
        day_in_year.append(time[n].timetuple().tm_yday)
    #    #day_in_year.append(time[n].timetuple().month)    
    day_in_year = np.array(day_in_year)
    
    if clean_masked :
        #def clean_masked() :
        #remove masked elements
        indexes_nonmasked = [o2.nonzero()]    
        o2 = o2[indexes_nonmasked].flatten()
        time =  time[indexes_nonmasked].flatten()
        day_in_year = day_in_year[indexes_nonmasked].flatten()        
        depth = depth[indexes_nonmasked].flatten()
    
    # take only depth < 2 m 
    indexes = np.where(depth < 5 )    
    o2 = o2[indexes].flatten()
    no3 = no3[indexes].flatten()
    po4 = po4[indexes].flatten()
    si = si[indexes].flatten()
    day_in_year = day_in_year[indexes].flatten()
    time = time[indexes].flatten()
    

    #arr = np.vstack((o2_sort2,day_in_year))
    #arr  = np.sort(arr, axis=1)
    if interpolate:
        #def to_interpolate() :
        #clean_masked()
        from scipy import interpolate
        f = interpolate.interp1d(day_in_year, o2, kind = 'nearest', # bounds_error=None, 
                             fill_value= 'extrapolate', assume_sorted=False)
    
        #f = interpolate.interp1d(time, o2_sort2, kind = 'nearest', # bounds_error=None, 
        #                         fill_value= 'extrapolate', assume_sorted=False)
        #from scipy.interpolate import CubicSpline
        #cs = CubicSpline(arr[1],arr[0])
        xnew = np.arange(1, 366, 1)
        ynew = f(xnew) 
    

        plt.plot( xnew, ynew, 'o--')
     
     
    fig = plt.figure()  
      
    ax = fig.add_subplot(221)
    ax1 = fig.add_subplot(222)   
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(224)       
     
    ax.set_title('O2')
    ax1.set_title('NO3')
    ax2.set_title('PO4')
    ax3.set_title('Si')    
    
    ax.plot(day_in_year,o2,'o')    
    ax1.plot(day_in_year,no3,'o') 
    ax2.plot(day_in_year,po4,'o')    
    ax3.plot(day_in_year,si,'o')     
    
    plt.show()

use_netcdf4()

'''
data = xr.open_dataset(file)
#data = data[:,0]
data.isnull() 
data = data.isel(N_SAMPLES = 0)

grouped_data = data.groupby('date_time.month')
#maxs = grouped_data.max() 
means = grouped_data.mean() 

#data.dims = ['var1','date_time']
#o2_time_series = np.array(data['var4'].sel(N_SAMPLES = 0)*44.6)
#print o2_time_series

#time_series = data['date_time.dayofyear']
#o2_time_series.values.tolist()


#print (o2_time_series.coords, '\nttt', time_series)
#o2_time_series.plot()
#plt.show()



#data = data.isel(N_SAMPLES = 0) 
#o2 = data['var4']*44.6 # convert to mM
#time = data['date_time.dayofyear']
#print (time)
#plt.plot(time,o2,'o')
#data.set_coords([ 'var1','date_time'])
#print (data.coords)
#fh = Dataset(file) 
#print (list(data.variables.keys()))



ds_min = grouped_data.min()
ds_max = grouped_data.max()
ds_mean = grouped_data.mean()

#x = ds_mean.var4.dayofyear

#print (x,ds_min.var4.dayofyear)
plt.plot(x, ds_max.var4*44.6 ,'ro')
plt.plot(x, ds_mean.var4*44.6 ,'mo')
plt.plot(x, ds_min.var4*44.6 ,'ko')

#ds_max.var4.plot()
#ds_mean.var4.plot()

plt.show()
#ds_min.var4.plot(kind='scatter',x ='date_time.month',y= ds_min)
#.reset_index()
'''

'''

#plt.scatter(ds_min.var4.keys, ds_min.var4.coords)
o2 = data['var4']
# first value 
#o2_surf = np.array((o2.isel(N_SAMPLES = 0)) * 44.6)
time = np.array(data['date_time'])
#time_surf = (time.isel(N_SAMPLES = 0))

#data.date_time.dt.month.plot() 
#plt.plot_date(time,o2_surf)
#o2_surf.plot.line('b^')
plt.show()

#print (o2_surf)
#time = data['date_time'] #.dt.month

#surface = data[(data['var1'] < 2)]     #(data['var1'] < 2).mean()
#print (data)

# take the shallowest part 

#grouped_data = data[data['date_time.month'] < 9]
#print (grouped_data)
#depth =  grouped_data['var1']
#mins = grouped_data.all() 
#o2min = mins.var4
#o2min = (o2min*44.6)

#maxs = grouped_data.max() 
#means = grouped_data.mean() 
#print(grouped_data.groups.keys())
#months = [3,4,5,7,8,9,11]

#plt.plot(months,o2min,'o')
#o2_surf.plot.line('b^')
#plt.show()
#count = grouped_data.count()
#times = num2date(time[:],time.units)





depth =  np.array((fh.variables['var1'][:]))
depth = ma.masked_invalid(depth)
depth = ma.masked_where(depth   == -10000000000.0 , depth )
o2 = np.array((fh.variables['var4'][:]))
o2 = ma.masked_invalid(o2)
o2 = ma.masked_where(o2  == -10000000000.0 , o2)
o2 = o2*44.6 # convert to mM
'''

'''
to_df = {'depth': [depth],'o2': [o2]} 
df = pd.DataFrame(to_df)
print (df.dtypes)

fig = plt.figure(figsize=(7,7))
for n in range(0,len(depth)):
    if depth[n,0] < 2: 
        plt.plot(o2[n,:2],depth[n,:2],'o--')

    
#plt.ylim(10,0)    
#plt.show()

'''   
        
        