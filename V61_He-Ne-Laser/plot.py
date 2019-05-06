import numpy as np
import scipy.optimize
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy import optimize
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy import stats
from matrix2latex import matrix2latex

def parabola(x, a, b, c):
    return a*x*x + b*x + c


#Gitter
#100 Spalte/mm -> 1*10^5 Spalte pro m
d_1=6.3 #cm
d_2=6.2 #cm
l=95.7 #cm

#HIER DIE AUSWERTUNG FÜR DIE WELLENLÄNGE




#plankonkav
print('plankonkav:')
l, I = np.genfromtxt('data/plankonkav.txt', unpack=True)

hr = ['$L$/cm', '$I$/µA']
m = np.zeros((21, 2))
m[:,0] = l
m[:,1] = I
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

plt.plot(l, I, 'rx', mew=0.5, label='Messwerte')
plt.xlabel(r'$L/$cm')
plt.ylabel(r'$I$/µA')
plt.xlim()
plt.ylim()
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/plankonkav.pdf')
plt.clf()

#plankonkav
print('konkavkonkav:')
l, I = np.genfromtxt('data/konkavkonkav.txt', unpack=True)

hr = ['$L$/cm', '$I$/µA']
m = np.zeros((29, 2))
m[:,0] = l
m[:,1] = I
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

plt.plot(l, I, 'rx', mew=0.5, label='Messwerte')
plt.xlabel(r'$L/$cm')
plt.ylabel(r'$I$/µA')
plt.xlim()
plt.ylim()
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/konkavkonkav.pdf')
plt.clf()

#TEM00
print('TEM00:')
s, I = np.genfromtxt('data/tem00.txt', unpack=True)

hr = ['$L$/cm', '$I$/nA']
m = np.zeros((63, 2))
m[:,0] = s
m[:,1] = I
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

plt.plot(s, I, 'rx', mew=0.5, label='Messwerte')
plt.xlabel(r'$s$/cm')
plt.ylabel(r'$I$/µA')
plt.xlim()
plt.ylim()
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/tem00.pdf')
plt.clf()

#TEM01
print('TEM01:')
s, I = np.genfromtxt('data/tem01.txt', unpack=True)

hr = ['$L$/cm', '$I$/nA']
m = np.zeros((63, 2))
m[:,0] = s
m[:,1] = I
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

plt.plot(s, I, 'rx', mew=0.5, label='Messwerte')
plt.xlabel(r'$s/$cm')
plt.ylabel(r'$I$/µA')
plt.xlim()
plt.ylim()
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/tem01.pdf')
plt.clf()

#polarisation
print('Polarisation:')
phi, I = np.genfromtxt('data/polarisation.txt', unpack=True)

hr = ['$L$/cm', '$I$/µA']
m = np.zeros((36, 2))
m[:,0] = phi
m[:,1] = I
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

plt.plot(phi, I, 'rx', mew=0.5, label='Messwerte')
plt.xlabel(r'$\phi/$°')
plt.ylabel(r'$I$/µA')
plt.xlim()
plt.ylim()
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/polarisation.pdf')
plt.clf()
