# Set of scripts for the paper 

### main.py 
Script to read all the nc files in a folder with different time periods from ROMS,
extracts data from one station, combines all time steps and writes to new netcdf. 
Calculates z in meters from sigmas 

### nc_to_each_day.py 
 script reads netcdf file created in main.py script, interpolates the  data to one day 
 time step. 
 Writes the data to new nectdf 
