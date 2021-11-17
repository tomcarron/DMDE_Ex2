'''
Python script for DMDE2 by Tom Carron.
'''
import numpy as np
import astropy
import astropy.constants as const
import astropy.units as u
import camb
from matplotlib import pyplot as plt
import os

'''
Problem 1: The matter power spectrum
a) Plot the matter power spectrum at the following epochs: z = 0.0, 0.5, 2.0, 10.0
Use CAMB and the parameters from the file planck 2018.ini .
Explain why the plot looks as it does.
'''
#The matter power spectrum describes the density contrast of the universe (the difference between the local density and the mean density) as a function of scale.
pars=camb.read_ini('planck_2018.ini')
results=camb.get_results(pars)
print("Calculating results...")
#print(pars)
epochs=np.array([0.0,0.5,2.0,10.0])
k=np.linspace(1e-4,10,1000)
pars.set_matter_power(epochs, kmax=10.0)
#pars.NonLinear=model.NonLinear_none
results=camb.get_results(pars)
kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
s8 = np.array(results.get_sigma8())

print("Plotting matter power spectra")
plt.figure(1)
for i, (epochs) in enumerate(epochs):
    plt.loglog(kh,pk[i,:],label='$z=$'+str(epochs))
plt.ylabel('$P(k,z)$ $[(h^{-1}Mpc)^3]$')
plt.xlabel('$k$ $[h Mpc^{-1}]$')
plt.legend()
plt.grid()















plt.show()


