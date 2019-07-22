import xarray as xr

roms_path = r"Data/ROMS_Laptev_Sea_NETCDF3_CLASSIC_east_var2.nc"
path_brom_in = r'C:\Users\elp\OneDrive - NIVA\BROM_linux_output\Baseline_B\water.nc'

ds = xr.open_dataset(roms_path)
depth = ds.depth.values
dzs = [round(abs(depth[k+1] - depth[k]),3) for k in range(len(depth)-1)]
print('ROMS ',dzs,depth) # dzs,'\n',depth)

ds_brom = xr.open_dataset(path_brom_in)
depth_brom = ds_brom.z.values
dzs_brom = [round(abs(depth_brom[k+1] - depth_brom[k]),3) for k in range(len(depth_brom)-1)]
print (' ')
print ('BROM',dzs_brom,depth_brom)


def make_df_sum(s):
    ''' Create scenario, sum rates of dissolution  from 
        different bubbles milliM/sec '''

    # Get rates of dissolution calculated in the Single Bubble Model
    df = pd.read_csv(bub_path,delimiter = '\t' )    
    sizes = s[0]
    df1 = df[df.radius == sizes[0]].reset_index()
    df_sum = df1.loc[:,['depth','rad_evol','met_cont','met_flow','vbub']]
    
    # Make sum according to scenario 
    for n,num in enumerate(sizes):          
        df2 = df[df.radius == num].reset_index()      
        for col in df_sum.columns[1:] :
            df_sum[col] = df_sum[col].add(df2[col],fill_value= 0)   
                      
    df_sum.met_cont = df_sum.met_cont*1000*0.5  # milliM
    df_sum.met_flow = df_sum.met_flow*1000*0.5  # milliM/sec  #50% of the time 
    return df_sum #Rates of dissolution for sum of bubbles for each scenario, milliM/sec