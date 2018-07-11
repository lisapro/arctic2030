import matplotlib.pyplot as plt

shakhova = [100,300]
black_sea = 800
svalbard = [100,300]
#ostenfjorden = 
w = 2
fig, ax = plt.subplots(figsize = (6,2))
ax.set_title(r'Bottom concentrations of dissolved CH$_4 $ nM')
ax.hlines(0,shakhova[0],shakhova[1],linewidth = w)
ax.annotate('ESAS (Shakhova et.al 2010)',(shakhova[0],0.2))
ax.hlines(1,black_sea,black_sea+20,linewidth = w)
ax.annotate('Black Sea \n(Schubert et al. 2006)',
            (black_sea-150,1.2))
ax.hlines(2,svalbard[0],svalbard[1],linewidth = w)
ax.annotate('Svalbard (Jansson et.al 2017)',(svalbard[0],1.7))
#plt.scatter(black_sea,1)
ax.set_yticks([])
plt.show()