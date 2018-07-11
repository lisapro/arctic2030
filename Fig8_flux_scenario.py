import matplotlib.pyplot as plt
import netCDF4
import xarray as xr


def plot_flux_scenario(sc_name):
    
    fig, axes = plt.subplots(nrows = 2,figsize = (6,4))
    ds = xr.open_dataset('Data\Laptev_average_year_2year.nc')
    f = ds[sc_name].loc['1991-1':'1991-12'] 
    
    cm = plt.get_cmap('Reds')
    
    f.plot(x='time',y = 'depth',ax = axes[0],
           cmap = cm,vmin = 0,vmax = 0.001)
    f.plot(x='time',y = 'depth',ax = axes[1],
           vmin = 0,vmax = 0.00001,cmap = cm) 
    
    axes[0].set_xticklabels([])
    axes[1].set_xticklabels(['Feb','Apr','June','Aug','Oct','Dec'])
    
    axes[0].set_xlabel(' ')
    axes[1].set_xlabel(' ')
    
    axes[0].set_ylim(10,0)
    axes[1].set_ylim(80,10)
    
    plt.show()

#plt.savefig('Data/fluxes_experiment.png',transparent = True)

if __name__ == '__main__':
    #plot_flux_scenario('B2_10f')
    plot_flux_scenario('B2_30f')