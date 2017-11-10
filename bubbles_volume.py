
# coding: utf-8

# In[1]:


#объем шара = 4/3 pi R^3


# In[2]:


import numpy as np


# In[3]:


'''
========================================
Calculation of the INITIAL BUBBLE VOLUME 
========================================
'''
# *1e-09 # volume in mm^3 to volume in meters^3
# *1e-03 # volume in mm^3 to volume in cm^3
# не домножать ни на что в выражении return функции vol_bub - получим мм^3

def vol_bub(r):
    return (4/3*np.pi*r**3)*1e-03 # объем газа в пузыре, cм^3 
    


# In[4]:


'''
=============================================
Calculation of the INITIAL BUBBLE MOL CONTENT
based on previous function
=============================================
'''
# подсчет начального кол-ва вещества (моль) на основании имеющихся радиусов (см) 
# с использованием уравнения Менделеева-Клайперона (PV = nRT)

def mol_bub():
    Pres = 7. # давление газа, атм ! что мы тут задаем????! давление газа как? или давление пузырька? или окружающей среды??
    Temp = 274.15 # температура, °К, сейчас здесь: 1°С (берем температуру придонной воды?)
    Rgas = 82.0575 # универсальная газовая постоянная, (atm*cm^3)/mol*°K
    vol = vol_bub(r) # объем газа на основании данных о радиусе, см^3
    return ((Pres*vol)/(Rgas*Temp))


# In[ ]:


'''
=============================================
Calculation of the INITIAL BUBBLE COMPOSITION
=============================================
if we assume that there are not only methane 
in the bubbles
---------------------------------------------
'''

def bub_content():
    met_bub=(mol_bub())/100*95
    o2_bub=(mol_bub())/100*3
    co2_bub=vol=(mol_bub())/100*2
    return(met_bub,o2_bub,co2_bub)
#print(met_bub, o2_bub, co2_bub)


# In[ ]:


volume =[] # пустые массивы # объем метана (см^3)
gas_bub =[] # общее кол-во газа (моли) (100%) в пузырьке - можем задавать как метан, так и различные соотношения.
met_bub =[] # содержание метана (моль) в пузыре
co2_bub =[] # содержание углекислого газа (моль) в пузыре
o2_bub =[] # содержание кислорода (моль) в пузыре

radius = np.arange(0.25,8.25,0.25) # создание набора радиусов (r, мм) с определенным шагом для подстчета формулой
for r in radius:
    vol_bub(r)
    mol_bub()
    bub_content()
    
    volume.append(vol_bub(r))
    gas_bub.append(mol_bub())
    met_bub.append(met_bub)
    co2_bub.append(co2_bub)
    o2_bub.append(o2_bub)
    
    #type(a)
    #print(r,vol_bub(r),mol_bub())
out_array = np.column_stack((radius, volume, mol_bub, met_bub,o2_bub,co2_bub))
np.savetxt('mol_res.dat', out_array, delimiter=' ') 
    


# In[37]:




