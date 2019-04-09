import pandas as pd
n = ['event','date','lat','lon','elevation','depth','ch4water','ch4 sed','kCH4','-','-','-','-']
df = pd.read_csv(r'C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Pangaea_data\Bussmann_2016.tab',skiprows = 47,delimiter = '\t',names = n)
#d = pd.Series(df.kCH4.values,dtype = 'float')
print(df.kCH4)