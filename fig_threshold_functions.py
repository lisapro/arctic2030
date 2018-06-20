import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

s_omso_no3 = 5.0 # no source
s_omso_o2 = 10 #15-30 find link
s_omch_so4 = 15000 # no source 

# K_DON_ch4 = 
DON = 6 #(6-11.25 range from 
#https://www.soest.hawaii.edu/oceanography/courses/OCN623/Spring2013/Organics-2013-handouts-2-14.pdf 
k_don_ch4 = 0.00014 #1/day
 
conc_no3 = np.arange(0,10,0.1)
conc_o2 = np.arange(0,50,0.1)
conc_so4 = np.arange(14990,15010,0.01)


def thresh_higher(conc,coef):
    thresh = 0.5 + 0.5*np.tanh(conc-coef) 
    return thresh

def thresh_lower(conc,coef):
    thresh = 0.5 - 0.5*np.tanh(conc-coef) 
    return thresh

r_no3 = thresh_lower(conc_no3,s_omso_no3)
r_o2 = thresh_lower(conc_o2,s_omso_o2)
r_anti_o2 = thresh_higher(conc_o2,s_omso_o2)
#r_so4 = 1.- 0.5*(1.+np.tanh(conc_so4-s_omch_so4))
#ax.text(6,0.5,'1.- 0.5*(1.+np.tanh(conc_no3-s_omso_no3))')
#r_om = (r_no3+r_o2+r_so4)*DON*k_don_ch4
fig = plt.figure(figsize = (4,5)) 
gs = gridspec.GridSpec(2,1)
gs.update(top = 0.9,bottom = 0.1, hspace = 0.4,left = 0.15,right = 0.9)
ax,ax1 = plt.subplot(gs[0]),plt.subplot(gs[1]) #,plt.subplot(gs[2])
#ax2 = ax.twiny()  # instantiate a second axes that shares the same x-axis
#,(ax,ax1,ax2) = plt.subplots(3, sharey = True) 

ax.axvline(s_omso_no3, linestyle='--', c='k',label=r'kNO3')
ax1.axvline(s_omso_o2, linestyle='--', c='k',label=r'kO2')
#ax2.axvline(s_omch_so4, linestyle='--', c='k',label='s_so4')

ax.plot(conc_no3, r_no3,'-',label='s_no3') 
ax1.plot(conc_o2,r_o2,'-',label='threshold O$_2$ lower')
ax1.plot(conc_o2,r_anti_o2,'-',label= 'threshold O$_2$ higher')

#ax2.plot(conc_so4,r_so4,'-',label='s_omso_so4',c = 'r')

ax.set_title('nitrate')
ax1.set_title('oxygen ')
#ax2.set_xlabel('sulfate ')
#ax2.legend()
ax.set_ylim(-0.1,1.1)
ax1.legend()

ax.legend(loc = 'center left')

plt.show()