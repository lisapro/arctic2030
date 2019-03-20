import numpy as np

#(4/3*np.pi*r**3)*1e-03 # Volume of gas in the bubble, cм^3
# подсчет начального кол-ва вещества (моль) на основании имеющихся радиусов (см) 
# с использованием уравнения Менделеева-Клайперона (PV = nRT)

def vol_mol_bub(r,t = -1):
    Pres = 7.9             # Pressure atm
    Temp = t + 273.15      # Temprature °К
    Rgas = 82.0575         # Gas constant cm3 atm K−1 mol−1 
    r = r/10               # radius, cm
    vol = (4/3*np.pi*r**3) # Volume of bubble cm^3
    #pv = nrT
    #methane content in the bubble millimole
    gas_mol=(Pres*vol)/(Rgas*Temp)*1000  
    return gas_mol

r8 = vol_mol_bub(8) 
r6 = vol_mol_bub(6) 
r5 = vol_mol_bub(5)
r4 = vol_mol_bub(4)
r3 = vol_mol_bub(3)
r2 = vol_mol_bub(2)
r1 = vol_mol_bub(1)

B1 = r4+r3+r2+r1
B2 = r3+r1
B3 = r3+r6+r2
print ('B1 SWI flux',B1, 'Millimole/sec')
print ('B2 SWI flux',B2, 'Millimole/sec')
print ('B3 SWI flux big bubbles', B3, 'Millimole/sec')
print ('8mm',r8, 'Millimole/sec')
print (1.15/1000,' micromol/m2 sec measured ')