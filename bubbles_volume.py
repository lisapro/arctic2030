# coding: utf-8
#объем шара = 4/3 pi R^3
import numpy as np

#Calculation of the INITIAL BUBBLE VOLUME 
# *1e-09 # volume in mm^3 to volume in meters^3
# *1e-03 # volume in mm^3 to volume in cm^3
# не домножать ни на что в выражении return функции vol_bub - получим мм^3

def vol_bub(r):
    return (4/3*np.pi*r**3)*1e-03 # объем газа в пузыре, cм^3 

# подсчет начального кол-ва вещества (моль) на основании имеющихся радиусов (см) 
# с использованием уравнения Менделеева-Клайперона (PV = nRT)
def mol_bub():
    Pres = 8. 
    # давление газа, атм ! что мы тут задаем????!
    # давление газа как? или давление пузырька? или окружающей среды??
    Temp = 274.15 # температура, °К, сейчас здесь: 1°С (берем температуру придонной воды?)
    Rgas = 82.0575 # универсальная газовая постоянная, (atm*cm^3)/mol*°K
    vol = vol_bub(r) # объем газа на основании данных о радиусе, см^3
    #n = PV/RT
    return ((Pres*vol)/(Rgas*Temp))

def bub_content():
    met_bub=(mol_bub())/100*95
    o2_bub=(mol_bub())/100*3
    co2_bub=vol=(mol_bub())/100*2
    return(met_bub,o2_bub,co2_bub)
#print(met_bub, o2_bub, co2_bub)

volume =[] # пустые массивы # объем метана (см^3)
gas_bub =[] # общее кол-во газа (моли) (100%) в пузырьке - можем задавать как метан, так и различные соотношения.
met_bub =[] # содержание метана (моль) в пузыре
co2_bub =[] # содержание углекислого газа (моль) в пузыре
o2_bub =[] # содержание кислорода (моль) в пузыре
r_to_file = []
radius = np.arange(0.25,8.25,0.25)
# создание набора радиусов (r, мм) с определенным шагом для подстчета формулой

for r in radius:
    v = vol_bub(r)
    m  = mol_bub()
    m_bub = bub_content()
    
    r_to_file.append(r)
    volume.append(v)
    gas_bub.append(m)
    met_bub.append(m_bub[0])
    co2_bub.append(m_bub[1])
    o2_bub.append(m_bub[2])

    #print(r,vol_bub(r),mol_bub())
out_array = np.column_stack((r_to_file, volume, met_bub,o2_bub,co2_bub))
import os 
dir_name = 'Data'   
script_dir = os.path.abspath(os.path.dirname(__file__))
dir_to_save = os.path.abspath(os.path.join(script_dir,dir_name))    
if not os.path.isdir(dir_to_save):
    os.makedirs(dir_to_save)
    
np.savetxt('{}\mol_res.dat'.format(dir_to_save), out_array, delimiter=' ') 
    