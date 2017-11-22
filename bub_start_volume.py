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

def vol_mol_bub(r):
    Pres = 7.5 # давление газа, атм ! что мы тут задаем????! давление газа как? или давление пузырька? или окружающей среды??
    Temp = 274.15 # температура, °К, сейчас здесь: 1°С (берем температуру придонной воды?)
    Rgas = 82.0575 # универсальная газовая постоянная, (atm*cm^3)/mol*°K
    vol = (4/3*np.pi*r**3)*1e-03
    # объем газа в пузыре, cм^3 на основании данных о радиусе
    # *1e-09 # volume in mm^3 to volume in meters^3
    # *1e-03 # volume in mm^3 to volume in cm^3
    # не домножать ни на что в выражении return функции vol_bub - получим мм^3
    gas_mol=(Pres*vol)/(Rgas*Temp) #содержание газа в пузырьке, моли
    return (r,vol,gas_mol) # если так: ([r,vol,mol]) - вернет tuple

dir_name = 'Data'   
script_dir = os.path.abspath(os.path.dirname(__file__))
dir_to_save = os.path.abspath(os.path.join(script_dir,dir_name))    
if not os.path.isdir(dir_to_save):
    os.makedirs(dir_to_save)

# In[26]:
mol_bub_list =[]
radius = np.arange(0.25,8.25,0.25) # создание набора радиусов (r, мм) с определенным шагом для подстчета формулой
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


