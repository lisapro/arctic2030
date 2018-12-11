import matplotlib.pyplot as plt
import netCDF4
import xarray as xr
import numpy as np
import matplotlib.ticker as mtick 
import seaborn as sns
sns.set()
brom_path = 'Data\Laptev_baseline.nc'

#cm = plt.get_cmap('GnBu') #('gist_stern_r') # #('Reds')
cm = plt.get_cmap('Reds')
    
def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
    
def Fig9_plot_flux_scenario(sc_name,save = False):    
    
    fig, ax =  plt.subplots(figsize = (8,2))
    
    ds = xr.open_dataset(brom_path)
    f = ds[sc_name].loc['1992-1':'1992-12'] 

    #f = f.where(f > 0)
    f = f.where(f == 0)
    cm = plt.get_cmap('Greys') #'Reds')
    format = mtick.FuncFormatter(fmt)
    ax.set_facecolor('k')

    cs = f.plot(x='time', y = 'depth', ax = ax,
                add_colorbar = False,levels = 5, 
            cmap = cm,vmin = 0,vmax = 0.0003) 

    #plt.colorbar(cs,ax = ax,format = format) 
         
    ax.set_xticklabels(['Feb','Apr','June',
                             'Aug','Oct','Dec'])

    ax.set_ylim(80,0)     
    ax.set_xlabel(' ')
    ax.set_ylabel('depth, m')           
    ax.set_title(sc_name)
    plt.tight_layout()#pad, h_pad, w_pad, rect
    
    if save == True:
        plt.savefig('Data/Random_seep_day.png'.format(sc_name),
                    transparent = False)
    else:    
        plt.show()

def plot_flux_scenario_one(sc_name,save = False):    
    
    fig, axes = plt.subplots(nrows = 2,figsize = (8,5))
    
    ds = xr.open_dataset(brom_path)
    f = ds[sc_name].loc['1992-1':'1992-12'] 

    #f = f.where(f > 0)
    #f = f.where(f == 0)
    cm = plt.get_cmap('Greys') #'Reds')


    format = mtick.FuncFormatter(fmt)
    axes[0].set_facecolor('w')
    axes[1].set_facecolor('w')
    cs = f.plot(x='time', y = 'depth', ax = axes[0],
                add_colorbar = False,levels = 50, 
            cmap = cm,vmin = 0,vmax = 0.0003)    
    plt.colorbar(cs,ax = axes[0],format = format) 
    
    cs2 = f.plot(x='time', y = 'depth', ax = axes[1],
           levels = 20,add_colorbar = False,
           vmin = 0,vmax = 0.00003,cmap = cm)
         
    plt.colorbar(cs2,ax = axes[1],format = format)    

    m = ['Feb','Apr','June','Aug','Oct','Dec'] 
    axes[0].set_xticklabels(m)
    axes[1].set_xticklabels(m)
    axes[0].set_ylim(5,0)
    axes[1].set_ylim(80,0)

    for n in [0,1]:        
        axes[n].set_xlabel(' ')
        axes[n].set_ylabel('depth, m')   
          
    axes[0].set_title(sc_name)
    
    
    if save == True:
        plt.savefig('Data/fluxes_experiment_{}.png'.format(sc_name),
                    transparent = False)
    else:    
        plt.show()

def plot_flux_scenarios(sc_name,sc_name2,save = False):    
    
    fig, axes = plt.subplots(nrows = 2,figsize = (6,4), dpi=150)
    
    ds = xr.open_dataset(brom_path)
    f = ds[sc_name].loc['1992-1':'1992-12'] 
    f2 = ds[sc_name2].loc['1992-1':'1992-12']    
    

    format = mtick.FuncFormatter(fmt)
    
    cs = f.plot(x='time', y = 'depth', ax = axes[0],
                add_colorbar = False,levels = 15, 
            cmap = cm,vmin = 0,vmax = 0.0005)    
    plt.colorbar(cs,ax = axes[0],format = format) 
    
    cs2 = f2.plot(x='time', y = 'depth', ax = axes[1],
           levels = 15,add_colorbar = False,
           vmin = 0,vmax = 0.00005,cmap = cm)
         
    plt.colorbar(cs2,ax = axes[1],format = format)    
     
    axes[0].set_xticklabels([])
    axes[1].set_xticklabels(['Feb','Apr','June',
                             'Aug','Oct','Dec'])
    axes[0].set_ylim(80,0)
    axes[1].set_ylim(80,0)
    
    for n in [0,1]:        
        axes[n].set_xlabel(' ')
        axes[n].set_ylabel('depth, m')   
           
    axes[0].set_title(r'Scenario B1_50 $ flux\  mmol\ CH_4 \cdot m^2 sec^{-1}$')
    axes[1].set_title(r'Scenario B2_50 $ flux\  mmol\ CH_4 \cdot m^2 sec^{-1}$')
    
    if save == True:
        plt.savefig('Data/Fig9_fluxes_experiments.png',
                    transparent = False)
    else:    
        plt.show()
              
if __name__ == '__main__':

    Fig9_plot_flux_scenario('B2_50f',save = False)

    #plot_flux_scenarios('B1_50f','B2_50f',False)

    #plot_flux_scenarios('B1_50f','B2_50f',False
