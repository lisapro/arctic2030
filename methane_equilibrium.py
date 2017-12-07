'''
Created on 30. nov. 2017

@author: ELP
'''

import math

from netCDF4 import Dataset
import tkinter as tk # python3
from tkinter.filedialog import askopenfilename # python 3
#To show only the dialog without any other GUI elements
root = tk.Tk()
root.withdraw()
import os
def calc_methane_surf(temp,sal,fg):
    
    ###          ONLY FOR SURFACE       ### 
    # Function calculates equilibrium methane concentration 
    # as a function of temperature and salinity for given 
    # atmospheric molecular fraction (fg) 
    # Reference: Denis A. Wiesenburg and Norman L. Guinasso, Jr.
    # Journal of Chemical and Engineering Data, Vol. 24, No. 4, 1979 
    # Equilibrium Solubilities of Methane,in Water and Sea Water 
    # Used equation (7) "a more general form of the atmospheric
    # equilibrium solubility equation that Weiss" 
    # Constants from table VI for nmol/L, from Moist Air at 1 Atm Total Pressure    

    abs_temp = temp + 273.15 # Absoulte Temp Kelvins
    # coefficients for nmol/l
    a1 = -415.2807  #-412.1710 nL/L
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


def calculate_vapor(temp,sal):
    #
    #function [vapor_press_atm] = vpress(S,T)
    # Copied from https://web.uvic.ca/~rhamme/
    '''% vpress   Vapor pressure of sea water
    %=========================================================================
    % vpress Version 2.0 : 27 October 2012
    %          Author: Roberta C. Hamme (University of Victoria)
    %
    % USAGE:  vapor_press = vpress(S,T)
    %
    % DESCRIPTION:
    %    Vapor pressure of sea water
    %
    % INPUT:  (if S and T are not singular they must have same dimensions)
    %   S = salinity    [PSS-78]
    %   T = temperature [degree C]
    %
    % OUTPUT:
    %   vapor_press = vapor pressure of seawater  [atm] 
    % 
    % AUTHOR:  Roberta Hamme (rhamme@uvic.ca)
    %
    % REFERENCE:
    %   Guide to Best Practices for Ocean CO2 Measurements
    %   Dickson, A.G., C.L. Sabine, J.R. Christian (Eds.) 2007
    %   PICES Special Publication 3, 191pp.
    %   Chapter 5: Physical and thermodynamic data
    %       Based on: Wagner, W., A. Pruss (2002) The IAPWS formulation 1995 
    %       for the thermodynamic properties of ordinary water substance for 
    %       general and scientific use, J. Phs. Chem. Ref. Data, 31, 387-535.
    %       AND Millero, F.J. (1974) Seawater as a multicomponent electrolyte 
    %       solution, pp.3-80.  In: The Sea, Vol. 5, E.D. Goldberg Ed.
    %
    % DISCLAIMER:
    %    This software is provided "as is" without warranty of any kind.  
    %=========================================================================
    
    % CALLER: general purpose
    % CALLEE: none
    
    %----------------------
    % Check input parameters
    %----------------------
    if nargin ~=2
       error('vpress.m: Must pass 2 parameters')
    end %if
    
    % CHECK S,T dimensions and verify consistent
    [ms,ns] = size(S);
    [mt,nt] = size(T);
    % Check that T&S have the same shape or are singular
    if ((ms~=mt) || (ns~=nt)) && (ms+ns>2) && (mt+nt>2)
       error('vpress: S & T must have same dimensions or be singular')
    end %if
    
    %------
    % BEGIN
    %------
    '''
    T = temp
    S = sal
    #%Calculate temperature in Kelvin and modified temperature for Chebyshev polynomial
    temp_K = T+273.15
    temp_mod = 1-temp_K/647.096;
    
    #%Calculate value of Wagner polynomial
    Wagner = -7.85951783*temp_mod +1.84408259*temp_mod**1.5 -11.7866497 * temp_mod**3 + 22.6807411*temp_mod**3.5 -15.9618719*temp_mod**4 + 1.80122502*temp_mod**7.5
    
    #%Vapor pressure of pure water in kiloPascals and mm of Hg
    vapor_0sal_kPa = math.exp(Wagner * 647.096 / temp_K)* 22.064 * 1000
    
    #%Correct vapor pressure for salinity
    molality = 31.998 * S /(1e3-1.005*S)
    osmotic_coef = 0.90799 -0.08992*(0.5*molality) +0.18458*(0.5*molality)**2 -0.07395*(0.5*molality)**3 -0.00221*(0.5*molality)**4
    vapor_press_kPa = vapor_0sal_kPa* math.exp(-0.018 * osmotic_coef * molality)
    
    #%Convert to atm
    vapor_press_atm = vapor_press_kPa/101.32501;

    return vapor_press_atm


def calc_methane_depth(temp,sal,fg,depth):
    
    # Bunsen is defined as the volume of gas
    # (corrected to standard temperature and pressure) absorbed in
    # a unit volume of water at the measurement temperature
    # when the partial pressure of the gas is 760 mm. The calculated
    # solubilities were corrected for the effect of dissolved gas
    # on the volume of the aqueous phase by using a value of 37
    # cm3 for the partial molal volume of methane (6').
    #
    # Coefficients and bunsen equation are from the paper 
    # 'Solubility of Methane in Distilled Water and Seawater' 1972 
    # Sachio Yamamoto,' James B. Alcauskas, and Thomas E. Crozier2
    # 
    abs_temp = temp + 273.15 # Absoulte Temp Kelvins
    a1 = -67.1962  #-68.8862 
    a2 = 99.1624 # 101.4956
    a3 = 27.9015 #28.7314 
    b1 = -0.072909 #0.076146
    b2 = 0.041674 #0.043970 
    b3 = -0.0064603 #-0.0068672
      
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
    c = bunsen * (p_tot - ((h/100.) * p_vapor)) * fg
    #print (temp,sal,bunsen,c) # to check 
    return c
import matplotlib.pyplot as plt

def call_met_profile():
    # function read nc file, calculates methage saturation with 2 var of 
    # functions and plots it 
    #We get the whole path to file
    file = askopenfilename(initialdir= os.getcwd(),
    filetypes =(("NetCDF file", "*.nc"),
    ("All Files","*.*")),title = "Choose a file.")
            
    fh = Dataset(file)     
    depth =  fh.variables['depth'][:]
    temp =  fh.variables['temp'][10,:] #10 random day, to change 
    sal =  fh.variables['sal'][10,:]
    
    methane = []
    methane2 = []
    
    fg = 1.8*10**(-6) # 1.87987 ppb to check 
    for n in range(0,len(temp)):
        t = temp[n]  
        s = sal[n]
        d = depth[n]
        met = calc_methane_surf(t,s,fg)
        met2 = calc_methane_depth(t,s,fg,d) * 44.6 * 10**6  #44.6 times nanomoles)  
        methane.append(met)
        methane2.append(met2)
        
        
    fig = plt.figure(figsize = (10,5))
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(1, 4)
    gs.update(wspace=0.3)
    ax = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    ax3 = plt.subplot(gs[3])  
      
    ax.plot(temp,depth,'o--')
    ax.set_title('Temperature C')
    
    ax1.plot(sal,depth,'go--')
    ax1.set_title('Salinity psu')
    
    ax2.plot(methane,depth,'ro--')
    ax2.set_title('Methane nM \nmethod 1 ')

    ax3.plot(methane2,depth,'ro--')
    ax3.set_title('Methane nM \nmethod 2')
   
    for axis in (ax,ax1,ax2,ax3):
        axis.set_ylim(80,0)
        
    plt.show()

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
    value = calc_methane_depth(temp,sal,fg,0) * 44.6 * 10**6 # (44.6 times nanomoles)  
    print ('test 2 ',value)
    test_bunsen = 0.04409

def test3():
    temp = 0 #C
    depth = 100 
    sal = 34 
    fg = 1.81*10**(-6) 
    value = calc_methane_depth(temp,sal,fg,depth) * 44.6 * 10**6 # (44.6 times nanomoles)  
    print (value)
    
call_met_profile()
#test2()
#test()
#test3()