import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
plt.style.use('ggplot') 
path = r'C:\Users\elp\OneDrive\Python_workspace\arctic2030\Data\all_for_sbm_79_mm_tab.dat'

df = pd.read_csv(path,'\t')
df['rad_evol_m'] = df['rad_evol']*10**-3
df['volume'] = 4/3*np.pi*df['rad_evol_m']**3  
df['concentrations'] = 1000*df['met_cont']/df['volume'] 


for n in np.arange(2.5,8,0.5): #)[1,2,3,4,5,6,7,8]:
    dff = df[df['radius'] == n]
    m = int(min(dff.concentrations))
    plt.plot(dff.concentrations,dff.depth,label = str(str(n)+'mm, min conc: '+str(m)))

#df4 = df[df['radius'] == 4]
#df1 = df[df['radius'] == 1]
#plt.plot(df8.concentrations,df8.depth,label = '8')

#plt.plot(df4.concentrations,df4.depth,label = '4')
#plt.plot(df1.concentrations,df1.depth,label = '1')
plt.title('Concentrations in bubbles Micromole/l')
plt.xlim(0,380000)
plt.ylim(80,0)
plt.legend()
plt.show()
#print (df.head())
