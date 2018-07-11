import roms_nc_seasonal_mean as roms
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np 


def plot_average_values():
    
    fig  = plt.figure(figsize=(9,6), dpi=100 )    
    gs0 = gridspec.GridSpec(1, 2)
    gs0.update(wspace = 0.3, top = 0.94,bottom = 0.1, 
               right = 0.9,left = 0.07) 

    #gs = gridspec.GridSpec(2, 2,)
    gs = gridspec.GridSpecFromSubplotSpec(2, 2,hspace=0.35,wspace = 0.05, 
                    width_ratios=(17, 1), subplot_spec=gs0[0])
    
    gs1 = gridspec.GridSpecFromSubplotSpec(2, 2,hspace=0.35,wspace = 0.05, 
                    width_ratios=(17, 1), subplot_spec=gs0[1])

    ax = fig.add_subplot(gs[0]) 
    cbax = fig.add_subplot(gs[1])
    ax1 = fig.add_subplot(gs[2])
    cbax1 = fig.add_subplot(gs[3])
    
    ax2 = fig.add_subplot(gs1[0])
    ax3 = fig.add_subplot(gs1[2])
    cbax3 = fig.add_subplot(gs1[3])
    
    for axis in (ax,ax1,ax3):
        axis.set_ylabel('depth, m')
        axis.set_ylim(80,0)   
         
    temp = roms.get_averaged_value('temp')
    sal = roms.get_averaged_value('sal')
    kz = roms.get_averaged_value('Kz_s')
    sm_hice, hice = roms.get_averaged_value_1d('hice')


    depth = temp.columns
    depth2 = kz.columns
    time = temp.index
    x,y = np.meshgrid(time,depth)
    x2,y2 = np.meshgrid(time,depth2)
    ax.set_title(r'Temperature, C')
    x_text = 0.05
    y_text = 1.15
    labels = ('(a) ','(b) ','(c)','(d)')
    for i,axis in enumerate((ax,ax1,ax2,ax3)):
        axis.text(x_text, y_text, labels[i], transform=axis.transAxes,
                 fontsize=14, fontweight='bold', va='top', ha='right')
        axis.set_xticks([1,31,60,91,121,152,182,213,244,274,305,335])
        axis.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    cax0 = ax.pcolormesh(x,y,temp.T,cmap = 'RdBu_r',vmin = -1.7,vmax = 1.7)
    cb = fig.colorbar(cax0, cax=cbax)
    #fig.colorbar(cax1,ax = c)
    ax1.set_title('Salinity, psu')
    cax1 = ax1.pcolormesh(x,y,sal.T)
    cb = fig.colorbar(cax1, cax=cbax1)
     #plt.title('Mean ice thickness,m ')
    ax2.set_title('Ice thickness, m')
    ax3.set_title('Salinity vertical diffusion\n coefficient, $m ^2 \cdot \sec$')
    ax2.plot(hice.index.values,hice.values,c = 'y')
    ax2.plot(hice.index.values,sm_hice,c='k')
    cax3 = ax3.pcolormesh(x2,y2,kz.T,cmap = plt.get_cmap('Reds'))
    cb = fig.colorbar(cax3, cax=cbax3)
    ax3.set_xlabel('Time')
    ax1.set_xlabel('Time')
    #plt.show()
    plt.savefig('Data/mean_roms_tshicekz.png',transparent = True)
    
plot_average_values()    
    