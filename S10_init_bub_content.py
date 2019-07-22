import numpy as np

#(4/3*np.pi*r**3)*1e-03 # Volume of gas in the bubble, cм^3
# подсчет начального кол-ва вещества (моль) на основании имеющихся радиусов (см) 
# с использованием уравнения Менделеева-Клайперона (PV = nRT)

def vol_mol_bub(r,t = -1):
    Pres = 8          # Pressure atm
    Temp = t + 273.15      # Temprature °К
    Rgas = 82.0575         # Gas constant cm3 atm K−1 mol−1 
    r = r/10               # radius, cm
    vol = (4/3*np.pi*r**3) # Volume of bubble cm^3
    #pv = nrT
    #methane content in the bubble millimole
    gas_mol=(Pres*vol)/(Rgas*Temp)*1000  
    return gas_mol

def flux(r):
    flux = vol_mol_bub(r,t = -1) / 10 
    return flux #(millimole/m2)

print ('Base',np.around(flux(4)*4,decimals = 3))
print ('IF',np.around(flux(4)*8,decimals = 3))
print ('SB',np.around(flux(2)*32,decimals = 3))
print ('RF',np.around(flux(4)*2,decimals = 3))

print (1.15/1000,' micromol/m2 sec measured ')
print (vol_mol_bub(4,t=-1)*4)
print (vol_mol_bub(4,t=0)*4)
print (vol_mol_bub(4,t=0)*4)
print (vol_mol_bub(4)*3) #Flux scenario
print (vol_mol_bub(4)*2) #Flux2 scenario
print (vol_mol_bub(4)*4) #Flux2 scenario