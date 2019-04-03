# coding: utf-8
# author: zagoven

import numpy as np
import os

# In[25]:


'''
=============================================
INITIAL BUBBLE VOLUME and MOL CONTENT,
COMPOSITION OF GASES
=============================================
'''
#(4/3*np.pi*r**3)*1e-03 # объем газа в пузыре, cм^3
# подсчет начального кол-ва вещества (моль) на основании имеющихся радиусов (см) 
# с использованием уравнения Менделеева-Клайперона (PV = nRT)

def vol_mol_bub(r,t = 1):
    Pres = 8 #7.5 #A Pressure atm
    Temp = t + 273.15 # Temprature °К
    Rgas = 82.0575 # Gas constant cm3 atm K−1 mol−1 
    r = r/10 # to cm
    vol = (4/3*np.pi*r**3) # Volume of bubble cм^3
    
    # объем газа в пузыре, cм^3 на основании данных о радиусе
    # *1e-09 # volume in mm^3 to volume in meters^3
    # *1e-03 # volume in mm^3 to volume in cm^3
    # не домножать ни на что в выражении return функции vol_bub - получим мм^3
    #pv = nrT
    #n = pv / rT
    gas_mol=(Pres*vol)/(Rgas*Temp) #содержание газа в пузырьке, моли
    
    #return (r,vol,gas_mol*1000) # если так: ([r,vol,mol]) - вернет tuple
    return (gas_mol*1000)

print ('B, Basic seep',np.around(vol_mol_bub(4,t = 1)*4,decimals=3)) 
print ('F, Reduced Flux ',np.around(vol_mol_bub(4,t = 1)*2,decimals=3)) 
print ('F, small Bubbles ',np.around(vol_mol_bub(2,t = 1)*32,decimals=3))
#print ((vol_mol_bub(4)+vol_mol_bub(3)+vol_mol_bub(2)+vol_mol_bub(1))*3.6/4)
def calculate_flux(r,n_bub):
    vol = (4./3.)*np.pi*r**3.*1e-06 #Liter  radius in mm liter   to liter **1e-06
    met_cont_sec = n_bub*vol /(22.413) # mole/sec
    met_cont_day = met_cont_sec*60*60*24  # Mole/sec
    #met_cont = [met_cont_sec,met_cont_day]
    return met_cont_sec 

def calculate_nbub(r,flux):
    vol = (4./3.)*np.pi*r**3.*1e-06 #liter
    n_bub = flux*22.4/vol
    return (np.around(n_bub,decimals=2))

#print ('ttt',vol_mol_bub(0.4)) 


#print (vol_mol_bub(4))
n_bub_init = 3.6 #bub
rad = 2 #mm
flux = calculate_flux(rad,n_bub_init)
#print (np.around(flux*60*60*24,decimals = 3) ,'Flux Mole CH4 /day;\n',np.around(flux*1.e3,decimals = 4),'Millimole /sec ' )
n_bub_calc = calculate_nbub(rad,flux)
#print (n_bub_calc,n_bub_calc == n_bub_init)

dir_name = 'Data'   
script_dir = os.path.abspath(os.path.dirname(__file__))
dir_to_save = os.path.abspath(os.path.join(script_dir,dir_name))    
if not os.path.isdir(dir_to_save):
    os.makedirs(dir_to_save)

# In[26]:
mol_bub_list =[]
radius = np.arange(0.25,8.25,0.25) 
# создание набора радиусов (r, мм) с определенным шагом для подстчета формулой



'''
for r in radius:
    vol_mol_bub = vol_mol_bub(r)
    mol_bub_list.append(vol_mol_bub)
mol_bub_array=np.array(mol_bub_list)
np.savetxt('{}\mol_vol_res.dat'.format(dir_to_save),mol_bub_array, delimiter=' ')
# In[31]:

#Если считаем все газовые составляющие(метан, углекислый, кислород):
met=95 #задаем процентное соотношение газа в пузырьке
co2=4
o2=1
met_bub=mol_bub_array[:,2]/100*met #обращение к 3 столбцу
co2_bub=mol_bub_array[:,2]/100*co2
o2_bub = mol_bub_array[:,2]/100*o2
results = np.column_stack((mol_bub_array,met_bub,co2_bub,o2_bub))


np.savetxt('{}\start_data.dat'.format(dir_to_save),results,delimiter=' ')

'''
