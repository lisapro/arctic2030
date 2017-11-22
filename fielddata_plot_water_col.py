'''
Created on 26. sep. 2017

@author: ELP
'''

from netCDF4 import Dataset
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy.ma as ma
import os,sys 
plt.style.use('ggplot')

cwd = os.getcwd()
# AMK 2015 September 

path_to_file = os.path.join(cwd,'Data\seep_data_amk.txt')
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
    o2 = r1[9]
    
# create figure with size close to a4 (vertical)
figure = plt.figure(figsize=(6, 8), dpi=100)
figure.suptitle('AMK 15 September, 2015', size = 13)
gs = gridspec.GridSpec(3, 2#,
                   #width_ratios=[1,1],
                   #height_ratios=[1,1]
                   )
gs.update(wspace=0.2,hspace = 0.3,left=0.1,
   right=0.97,bottom = 0.05, top = 0.9) 

color = '#74a4c7'
mcol = '#11202b'
ax00 = figure.add_subplot(gs[0])
ax01 = figure.add_subplot(gs[1])    
ax02 = figure.add_subplot(gs[2])
ax03 = figure.add_subplot(gs[3]) 
ax04 = figure.add_subplot(gs[4])
ax05 = figure.add_subplot(gs[5])    

ax00.set_title('Temperature')
ax00.plot(temp1,depth1,'o-', color = color,markerfacecolor = mcol)
ax00.plot(temp2,depth2,'o-',color = color,markerfacecolor = mcol)
ax00.set_ylim(70,0)


ax01.set_title('salinity')
ax01.plot(sal1,depth1,'o-',c= color,markerfacecolor = mcol)
ax01.plot(sal2,depth2,'o-',c= color,markerfacecolor = mcol)


ax02.set_title('NO3')
ax02.plot(no3_1,depth1,'o-',c= color,markerfacecolor = mcol)
ax02.plot(no3_2,depth2,'o-',c= color,markerfacecolor = mcol)


ax03.set_title('Si')
ax03.plot(si_1,depth1,'o-',c= color,markerfacecolor = mcol)
ax03.plot(si_2,depth2,'o-',c= color,markerfacecolor = mcol)


ax04.set_title('PO4')
ax04.plot(po4_1,depth1,'o-',c= color,markerfacecolor = mcol)
ax04.plot(po4_2,depth2,'o-',c= color,markerfacecolor = mcol)

ax05.set_title('alkalinity')
ax05.plot(alk_1,depth1,'o-',c= color,markerfacecolor = mcol)
ax05.plot(alk_2,depth2,'o-',c= color,markerfacecolor = mcol)

for axis in (ax00,ax01,ax02,ax03,ax04,ax05):
    axis.set_ylim(70,0)
   

plt.show()  