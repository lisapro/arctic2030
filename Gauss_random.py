# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 15:20:38 2017

@author: zagoven
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import random

#Draw samples from the distribution:
#mu, sigma = 4., 1.87 # mean and standard deviation P.S. (st deviation^2=variance)

radius = np.arange(1.,7.5,0.5)
mu = np.mean(radius)
sigma=np.std(radius)
s = np.random.normal(mu, sigma, 1000)

#np.savetxt('s.dat',s,delimiter=' ')
#Verify the mean and the variance:
abs(mu - np.mean(s)) < np.sqrt(sigma)
abs(sigma - np.std(s, ddof=1)) < np.sqrt(sigma)

count, bins, ignored = plt.hist(s,11, normed=True)
norm_function=1/(sigma * np.sqrt(2 * np.pi )) * np.exp( - (bins - mu)**2 / (2 * sigma**2))
plt.plot(bins,norm_function,linewidth=2, color='r')
plt.xlim(0.25,7.45)
#plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * sigma**2) ),linewidth=2, color='r')
plt.savefig('gauss.png', dpi=400, facecolor='w', edgecolor='r',
        orientation='portrait', papertype=None, forma=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)
plt.show()
