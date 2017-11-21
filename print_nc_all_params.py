
from netCDF4 import Dataset

fname = 'ROMS_Laptev_Sea_NETCDF3_CLASSIC_south_each_day.nc'
fh = Dataset(fname)   
names_vars = [] 
for names,vars in fh.variables.items():
    names_vars.append(names)  
    
print (names_vars)   