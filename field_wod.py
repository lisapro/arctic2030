'''
Created on 2. okt. 2017
#!/usr/bin/python
# -*- coding: utf-8 -*-
@author: ELP
'''
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

'''
names_vars = [] 
for names,vars in fh.variables.items():
    names_vars.append(names)  

    try:
        x = fh.variables[str(names)].long_name
        print(names, x) 
    except AttributeError :     
        print ('error')
'''        

fh = Dataset(r'C:\Users\ELP\workspace\arctic2030\src\data_from_WOD_COLLECTION_(2017-08-21T16-04-41).nc')
fh_roms = Dataset(r'ROMS_Laptev_Sea_NETCDF3_CLASSIC_south_each_day.nc')     
fh_water_brom = Dataset(r'water.nc')     
#d = {'depth': [depth],'o2': [o2],'time': [jd]} 

#datatypes 
#, 'col2': [3, 4]
#df2 = pd.DataFrame(data = d)
#df2[['o2']] = df2[['o2']].apply(pd.to_numeric)
#df = pd.DataFrame(df2, dtype='float')
#print (df2.dtypes)
#df2.plot(x = 'time',y ='depth') 
#hs = pd.Series(depth[:,0],index=jd)

fig = plt.figure(figsize=(7,7))

#hs.plot('o') #,title='%s at %s' % (o2.long_name,fh.id))


depth =  (fh.variables['var1'][:])
o2 = (fh.variables['var4'][:][:])
o2 = o2*44.6 # convert to mM
no3 = fh.variables['var7'][:][:]
po4 = fh.variables['var5'][:][:]
si = fh.variables['var6'][:][:]
alk = fh.variables['var12'][:][:]
alk = alk * 1000  
depth_roms = fh_roms.variables['depth'][:]
time_roms = fh_roms.variables['time'][:] 
depth_roms_2d = []
for n in enumerate(time_roms):
    depth_roms_2d.append(depth_roms) 
      
o2_roms = fh_roms.variables['o2'][:]     
no3_roms = fh_roms.variables['no3'][:][:]    
po4_roms = fh_roms.variables['po4'][:][:]   
si_roms = fh_roms.variables['Si'][:][:]
alk_roms = fh_roms.variables['Alk'][:][:] 
 
    
depth_water_brom = fh_water_brom.variables['z'][:]
o2_water_brom = fh_water_brom.variables['B_BIO_O2'][:]  
time_water_brom = fh_water_brom.variables['time'][:] 
depth_water_brom_2d = []
for n in enumerate(time_water_brom):
    depth_water_brom_2d.append(depth_water_brom)
    
no3_water_brom = fh_water_brom.variables['B_NUT_NO3'][:]  
po4_water_brom = fh_water_brom.variables['B_NUT_PO4'][:]  
si_water_brom = fh_water_brom.variables['B_NUT_Si'][:]  
alk_water_brom = fh_water_brom.variables['B_C_Alk'][:]
   
time =  fh.variables['date_time']
times = num2date(time[:],time.units)
sal = fh.variables['var3'][:][:]
temp = fh.variables['var2'][:][:]



no2 = fh.variables['var8'][:][:]
# pH = fh.variables['var9'][:][:] no pH in roms
# chl = fh.variables['var10'][:][:] no data 
plankt_biomass = fh.variables['var11'][:][:]

Ntot = fh.variables['var13'][:][:] #no2+no3
pCO2 = fh.variables['var14'][:][:] 


import os,csv
#AMK 2015 
cwd = os.getcwd()
# AMK 2015 September 


'''
data = pd.read_csv('results\seep_data_amk.txt', sep="\t")
x = data.depth < 5 
print (data[x]) #,data.O2ml)
data.style.highlight_null(null_color='red')
#pd.set_option('display.html')
print (data.columns)
#plt.scatter(data)
plt.ylim(100,0)
#plt.show()
'''


path_to_file = os.path.join(cwd,'results\seep_data_amk.txt')
with open(path_to_file, 'r') as f:
    # important to specify delimiter right 
    reader = csv.reader(f, delimiter='\t')
    r = []
    for row in reader:
        r.append(row)        
    r1 = np.transpose(np.array(r[1:]) ) # skip  header line
    
    depth1 = np.array(r1[6][0:6]).astype(float) 
    depth1 = depth1+71
    depth2 = np.array(r1[6][6:-6]).astype(float) 
    depth2 = depth2+71
    temp1 = r1[7][0:6]
    temp2 = r1[7][6:-6]
    sal1 = r1[8][0:6]
    sal2 = r1[8][6:-6]

    alk_1 = np.array(r1[11][0:6]).astype(float) 
    alk_1 = alk_1*1000.
    alk_2 = np.array(r1[11][6:-6]).astype(float)     
    alk_2 = alk_2 * 1000  
    po4_1 = r1[12][0:6]
    po4_2 = r1[12][6:-6]      
    si_1 = r1[13][0:6]
    si_2 = r1[13][6:-6]     
    no3_1 = r1[14][0:6]
    no3_2 = r1[14][6:-6]   

def plot_var(var,var2,var3,amkvar1,amkvar2,title):
    plt.clf()
    ax = fig.add_subplot(111)


    plt.scatter(var3,depth_water_brom_2d,alpha = 1,
               label= 'ipbm_brom',s = 3, c = "#a5a7aa") #  edgecolor = '#1c254a',
    plt.scatter(var2,depth_roms_2d,alpha = 0.7, #c = "w",
                label = 'roms',s = 3 )  #edgecolor = '#5a2d18',    
    plt.scatter(var,depth, label = 'wod',s = 4,c = '#a00a3d')   #edgecolor = 'k',    
    try:   
        plt.scatter(amkvar1,depth1,alpha = 1, c = 'k',
               label= 'AMK 2015',s = 5) #  edgecolor = '#1c254a',   
        plt.scatter(amkvar2,depth2,alpha = 1, c = 'k',s = 5) #  edgecolor = '#1c254a',         
    except ValueError: 
        pass 
                     
    ax.set_ylabel('depth (m)')
    plt.title(title)                
    plt.ylim(100,0)

    plt.legend()
    plt.show()
    #plt.savefig('results/'+str(title)+'.png')   
    

#plot_var(o2,o2_roms,o2_water_brom,None,None,'o2')    
#plot_var(no3,no3_roms,no3_water_brom,None,None,'no3')  
plot_var(po4,po4_roms,po4_water_brom,po4_1,po4_1,'po4')  
#plot_var(si,si_roms,si_water_brom,None,None,'si') 
#plot_var(alk,alk_roms,alk_water_brom,alk_1,alk_2,'Alk') '''
