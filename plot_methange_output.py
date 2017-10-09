'''
Created on 9. okt. 2017

@author: ELP
'''

from netCDF4 import Dataset
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


file = r'E:\Users\ELP\Fortran\bubl\out.dat '
with open(file , 'r') as f:
    ## important to specify delimiter right 
    #reader = csv.reader(f, delimiter=',')
    r = []
    for row in f:
        # if you don't know which delimiter is used 
        # print one row to view it 
        row = row.split()
        #print (row)
        #break
        r.append(row)        
    r1 = np.transpose(np.array(r[:])) 
    # time,zr/100.,rn/1000.,mn, kb
time = r1[0]
depth_bubble = r1[1] # z of bubble (m)?? 
radius_bubble = np.array(r1[2]).astype(float)  # r of bubble ! r in micrometer ?? 
s = np.array(radius_bubble - min(radius_bubble))*1000.
mn_volume  = np.array(r1[3]).astype(float)   # volume of methane 
kb = r1[4] # k mass transfer ! cm/sec
    
ymax = max(depth_bubble)
ch4_max = max(mn_volume)
ch4_min = min(mn_volume)

def plot_fig():
    figure = plt.figure(figsize=(8,6), dpi=100)
    gs = gridspec.GridSpec(1, 2)
    gs.update(wspace=0.2,hspace = 0.3,left=0.1,
       right=0.97,bottom = 0.05, top = 0.85)
    ax00 = figure.add_subplot(gs[0])
    ax01 = figure.add_subplot(gs[1])   
    for n in range(0,len(radius_bubble),100) :
        ax00.scatter(radius_bubble[n],depth_bubble[n],s = s[n],alpha = 0.2,c = 'k')
        ax01.scatter(mn_volume[n],depth_bubble[n],alpha = 0.2,c = 'r',s = s[n])        
    ax00.set_ylim(100,0)
    ax01.set_ylim(100,0)
    ax01.set_xlim(ch4_min,ch4_max)
    #ax01.plot(var2,depth)   
     
    plt.show()

    
plot_fig()    
