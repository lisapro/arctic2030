'''
Created on 30. nov. 2017

@author: ELP
'''

import math
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import tkinter as tk 
from tkinter.filedialog import askopenfilename 
root = tk.Tk()
root.withdraw()
import os
import matplotlib.gridspec as gridspec
    
def calc_methane_surf(temp,sal,fg):
    
    '''        ONLY FOR SURFACE      
    Function calculates equilibrium methane concentration 
    as a function of temperature and salinity for given 
    atmospheric molecular fraction (fg) 
    Reference: Denis A. Wiesenburg and Norman L. Guinasso, Jr.
    Journal of Chemical and Engineering Data, Vol. 24, No. 4, 1979 
    Equilibrium Solubilities of Methane,in Water and Sea Water 
    Used equation (7) "a more general form of the atmospheric
    equilibrium solubility equation that Weiss" 
    Constants from table VI for nmol/L, from Moist Air at 1 Atm Total Pressure 
    '''   

    abs_temp = temp + 273.15 # Absoulte Temp Kelvins
    a1 = -415.2807  # coefficients for nmol/l 
    a2 = 596.8104
    a3 = 379.2599
    a4 = -62.0757 
    b1 = -0.059160
    b2 = 0.032174 
    b3 = -0.0048198     
    
    ln_C = math.log(fg) + a1 + a2 * (100. / abs_temp) + a3 * math.log(
    abs_temp / 100.)  + a4 * (abs_temp / 100.) + sal * (
        b1 + b2*(abs_temp / 100.) + b3 * ((abs_temp / 100)**2))
    
    c = math.e ** ln_C #nmol/l 
    
    return c

def calculate_vapor(temp,S):

    '''Vapor pressure of sea water
    Matlab code Copied from https://web.uvic.ca/~rhamme/
    Author: Roberta C. Hamme 
    (University of Victoria,rhamme@uvic.ca)
    and changed from matlab to python 
        
    DESCRIPTION:  Vapor pressure of sea water
    S = salinity    [PSS-78]
    T = temperature [degree C]
    OUTPUT:vapor_press = vapor pressure of seawater  [atm] 
    REFERENCE: Guide to Best Practices for Ocean CO2 Measurements
    Dickson, A.G., C.L. Sabine, J.R. Christian (Eds.) 2007
    PICES Special Publication 3, 191pp.
    '''
    
    #Calculate temperature in Kelvin and modified 
    # temperature for Chebyshev polynomial
    temp_K = temp + 273.15
    temp_mod = 1-temp_K/647.096;
    
    #Calculate value of Wagner polynomial
    Wagner = (-7.85951783*temp_mod +1.84408259*temp_mod**1.5 - 
              11.7866497 * temp_mod**3 + 
              22.6807411*temp_mod**3.5 -15.9618719*temp_mod**4 +
              1.80122502*temp_mod**7.5)
    
    #Vapor pressure of pure water in kiloPascals and mm of Hg
    vapor_0sal_kPa = math.exp(Wagner * 647.096 / temp_K)* 22.064 * 1000
    
    #Correct vapor pressure for salinity
    molality = 31.998 * S /(1e3-1.005*S)
    osmotic_coef = (0.90799 -0.08992*(0.5*molality) +
                    0.18458*(0.5*molality)**2 -
                    0.07395*(0.5*molality)**3 -
                    0.00221*(0.5*molality)**4)
    vapor_press_kPa = vapor_0sal_kPa * math.exp(
                      -0.018 * osmotic_coef * molality)
    
    #Convert to atm    
    vapor_press_atm = vapor_press_kPa/101.32501 
    return vapor_press_atm


def calc_methane_depth(temp,sal,fg,depth):
    
    ''' Calculate equilibrium methane concentration
        Bunsen  -  the volume of gas (corrected to st.temp 
        and pressure) absorbed in a unit volume of water 
        at the measured temp when the part.pressure of the gas is 760 mm.
        
        Coefficients and bunsen equation are taken from the paper 
       'Solubility of Methane in Distilled Water and Seawater' 1972 
       Sachio Yamamoto,' James B. Alcauskas, and Thomas E. Crozier
    '''
    
    abs_temp = temp + 273.15 # Absoulte Temp Kelvins
    a1 = -67.1962  
    a2 = 99.1624 
    a3 = 27.9015 
    b1 = -0.072909
    b2 = 0.041674  
    b3 = -0.0064603 
      
    ln_bunsen = a1 + a2*(100. / abs_temp) + a3 * math.log(abs_temp / 100.) + (
         sal * ( b1 + b2 * (abs_temp / 100.)  + b3 * ((abs_temp / 100.)**2)) )
    
    bunsen = math.e ** ln_bunsen
        
    # p_tot - the total pressure (atm) 
    # fg -  the mole fraction of gas (fG) in the dry atmosphere
    # h is the relative humidity (percent)
    # p_vapor is vapor pressure of the solution (atm).
    
    try:
        import seawater as sw
        st_lat = 76.47
        p_dbar = sw.pres(depth,st_lat) # dbars
        p_tot = 0.987 * p_dbar/10. # convert dbars to atm 
    except ModuleNotFoundError:
        p_tot = 1 + depth/10.3 # in atm
            
    h = 100
    p_vapor = calculate_vapor(temp,sal)
     
    # Equation 
    c = (bunsen * (p_tot - ((h/100.) * p_vapor)) * 
        fg * 44.6 * 10**6) #(nanomoles)  
    
    return c,bunsen
  
def call_met_profile():
    
    ''' reads nc file with temperature and 
        salinity profiles, calculates methane 
        returns depth,temp,sal, 
        methane equilibrium solubility '''
    
    file = r'Data\Laptev_average_year_3year.nc'            
    fh = Dataset(file)     
    
    # Here we take random day in a year 
    # Since ts does not change much and do not 
    # influence solubility significantly
    
    depth =  fh.variables['depth'][:]    
    temp =  fh.variables['temp'][10,:]  
    sal =  fh.variables['sal'][10,:]
    
    methane = []        
    fg = 1.8*10**(-6) # 1.87987 ppb to check 
    for n in range(0,len(temp)):
        t = temp[n]  
        s = sal[n]
        d = depth[n]
        met = calc_methane_depth(t,s,fg,d)[0] 
        methane.append(met)
     
    return depth,temp,sal,methane 

def plot_eq_methane(): 
    
    depth,temp,sal,methane = call_met_profile()    
    fig = plt.figure(figsize = (3,5))

    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.3)
    ax = plt.subplot(gs[0])
    
    ax.plot(methane,depth,'ro--')   
    ax.set_title('Methane solubility nM')  
    ax.set_ylim(80,0)
    ax.set_xlim(0,35)    
    plt.show()       

def plot_ts_eq_methane():
    
    depth,temp,sal,methane = call_met_profile()
    fig = plt.figure(figsize = (10,5))
    gs = gridspec.GridSpec(1, 3)
    gs.update(wspace=0.3)
    ax = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
     
    ax.plot(temp,depth,'o--')
    ax.set_title('Temperature C')    
    ax1.plot(sal,depth,'go--')
    ax1.set_title('Salinity psu')   
    ax2.plot(methane,depth,'ro--')
    ax2.set_title('Methane nM \nmethod 2')
   
    for axis in (ax,ax1,ax2):
        axis.set_ylim(80,0)
        
    plt.show()

def calculate_flux(windspeed,ch4_water,temp,sal,depth,pCH4_air):
    
    ''' Function calculates air-sea flux of CH4
    Flux in [muM m**-2 d**-1] or flux_sec in [muM m**-2 sec**-1] 
    based on Wanninkhof et al. (2009) and L. Lapham et al.(2017)
    input parameters: 
    windspeed at 10 m [m*sec**-1]
    pCH4_air [ppm] = [uatm]  
    ch4_water - CH4 concentration at surface water layer,[nM/l] 
    temp [C], sal [psu], depth [m] '''
    
    bunsen = calc_methane_depth(temp,sal,pCH4_air,depth)[1] # [ml CH4 / ml water]  
      
    #p = nRt/v = c[nM/l] * R[atm*l*(mol**-1)*(K**-1)] * T[Kelvin]
    pCH4_water = ch4_water* 0.082057 * (temp + 273.15) * (10**-3) *22.4    #![uatm]
        
    # Coefficients from  Wanninkhof (2014) Table 1
    Sc_ch4 = (2101.2 - 131.54 * temp + 4.4931 * temp **2 
              - 0.08676 * temp **3 + 0.00070663 * temp ** 4) 
    s = bunsen / 22.4  # mole/l                           

    #coef 0.24 includes convertions from m/s to cm/h 
    k = ( 0.24 * (windspeed**2 )) * ((Sc_ch4/ 660.)**-0.5) # [cm/h] 
    flux = k * s * (pCH4_water - pCH4_air) * 24 / 100 # [muM m**-2 d**-1]
    flux_sec = flux / 86400. # [muM m**-2 sec**-1]
    return flux  


'''
def test(): 
    temp = 10 #Celsius
    sal = 34 
    # Molecular Fraction of Atmospheric Gases Used for Test Value Calculations, ppm     
    fg = 1.41*10**(-6) 
    value = calc_methane_surf(temp,sal,fg)  
    test_value_nmol = 2.146
    print ('test1',test_value_nmol, value)
 
        
def test2(): 
    temp = 10 #0.73  #Celsius
    sal = 34 #33.629 # commented values for bunsen test
    # Molecular Fraction of Atmospheric Gases Used for Test Value Calculations, ppm     
    fg = 1.41*10**(-6) 
    value = calc_methane_depth(temp,sal,fg,0) 
    conc = value[0]
    bunsen = value[1]
    print ('test 2 ',value)
    test_bunsen = 0.04409

def test3():
    temp = 0 #C
    depth = 100 
    sal = 34 
    fg = 1.81*10**(-6) 
    value = calc_methane_depth(temp,sal,fg,depth)
    conc = value[0] 
    bunsen = value[1]
    print (value)
  
  
def convert():
    # 22.4 liter - 1 mole 
    # 1 liter - x 
    # x = ( (1 [mole] * input[liter]) / 22.4 liter )[mole]
    input = 3. *10 **(-5) # ml/l 
    input = input * 10**(-3) # convert ml to l
    x  = ( 1 * input ) / 22.4 #[mole/l]
    x = x * 10 ** 9 #[nanomole/l]
    print (x, 'nM/l')   
    return x '''
 
     
#convert()


if __name__ == '__main__':

    #plot_ts_eq_methane()
    
    # test values from table 1  L. Lapham et al.(2017)  
    '''test_1line = calculate_flux(windspeed = 4.53, ch4_water = 13.69,
                     temp = 1.31,sal = 30.18, depth = 5.05,pCH4_air = 1.89)     
    print (test_1line)
    test_2line = calculate_flux(windspeed = 5.52, ch4_water =22.73 ,
                   temp = -0.66,sal = 29.11, depth = 5.03,pCH4_air = 1.89) '''
    