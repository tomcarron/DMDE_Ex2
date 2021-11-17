'''
Python script for DMDE2 by Tom Carron.
'''
import numpy as np
import astropy
import astropy.constants as const
import astropy.units as u
import astropy.cosmology as cosmo
from astropy.cosmology import Planck18_arXiv_v2 as cosmo
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
k=np.linspace(1e-4,1,1000)
k_h=k/cosmo.h
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

'''
Plot the current epoch matter power spectrum P(k,z=0) using the approximation given by equations 2.109-2.112 in the lecture notes (Bardeen et al. 1986). You may
use the present day cosmological parameters you found in the first exercise sheet from
Planck 2018.
'''

#2.112
shape_parameter=((cosmo.Om0)*(cosmo.h)*((2.7/cosmo.Tcmb0)**2)*(np.exp(-cosmo.Ob0-np.sqrt((cosmo.h)/0.5)*(cosmo.Ob0/cosmo.Om0)))).value
#2.111

def q(k):
    q = k / (shape_parameter*cosmo.h) #k in Mpc
    return q

#2.110
def transfer_function(k):
    T=(np.log(1+2.34*q(k))/2.34*q(k)) * ((1+(3.89*q(k))+((16.1*q(k))**2)+((5.49*q(k))**3)+((6.71*q(k))**4))**(-0.25))
    return T

#2.109
def P(k):
    n=1  #primordial power spectral index, n, has a value close to 1 (e.g., Komatsu et al. 2008)
    A=1 # The normalization can only be determined through observations
    P=A*(k**n)*(transfer_function(k))**2
    return P


#print(cosmo)

plt.figure(2)
plt.loglog(k_h,P(k_h),label='z=0')
plt.title('Power spectrum')
plt.xlabel('$k$')
plt.ylabel('$P(k,z=0)$')

#print(P(0))
'''
Probably a problem with units, too small and missing half the function
'''
'''
ok so i now think the problem is with the firsdt plot not the second one. See Schenider book and Cell 32 in the Jupyter notebook on this
'''


plt.show()


