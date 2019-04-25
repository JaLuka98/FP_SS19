import numpy as np
import scipy.optimize
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy import optimize
import matplotlib.pyplot as plt
from uncertainties import ufloat
from uncertainties import correlated_values
from matrix2latex import matrix2latex

#Zylinderkette
def linfit(x,a,b):
    return a*x+b

zylinder, f1, a1 , b1, A1, f2, a2, b2, A2 = np.genfromtxt('data/roehre.txt', unpack=True)

deltaf=np.log(1000*(f2-f1))
zylinderanzahl=np.log(zylinder)
linspace = np.linspace(-0.1, 2.7, 500)
params, covariance_matrix = optimize.curve_fit(linfit, zylinderanzahl, deltaf)
a, b = correlated_values(params, covariance_matrix)
print('Fit zur Schallgeschwindigkeitsbestimmung')
print('a=', a)
print('b=', b)

plt.plot(zylinderanzahl, deltaf, 'rx', mew=0.5, label='Messwerte')
plt.plot(linspace, linfit(linspace, *params),'b-', label='Ausgleichrechnung', linewidth=0.5)
#plt.yscale('log')
#plt.xscale('log')
plt.xlabel(r'$\log$(Anzahl der Zylinder)')
plt.ylabel(r'$\log(\Delta f)$')
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/roehre.pdf')
plt.clf()

zylinder=np.loadtxt('data/zyldat.txt', dtype=np.str)

for i in range (0, len(zylinder)):
    f, A = np.genfromtxt(zylinder[i], unpack=True)
    #plt.plot(f, A, 'ro', mew=0.5)
    plt.plot(f, A, 'b-', linewidth=0.5, label='Gemessenes Frequenzspektrum')
    f_osz=np.genfromtxt('data/oszispektrum.txt', usecols=(i))
    for a in range(0, len(f_osz)):
        plt.axvline(f_osz[a]*1000, color='r', linestyle='--', linewidth=0.4)
    plt.xlabel(r'$f/$Hz')
    plt.ylabel(r'$A$')
    plt.xlim(0, np.max(f))
    plt.ylim()
    plt.tight_layout()
    plt.legend()
    plt.grid()
    plt.savefig('build/zyl'+str(i+1)+'.pdf')
    plt.clf()


#H-Atom

f, A=np.genfromtxt('data/hatom/hatom180alles.dat', unpack=True)
f_osz=np.genfromtxt('data/hatom.txt', unpack=True)
plt.plot()
plt.plot(f, A, 'b-', linewidth=0.5, label='PC')
for i in range(0, len(f_osz)):
    plt.axvline(f_osz[i]*1000, color='r', linestyle='--', linewidth=0.5)
plt.xlabel(r'$f/$Hz')
plt.ylabel(r'$A$')
plt.xlim()
plt.ylim()
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/hatomalles.pdf')
plt.clf()



wasserstoff=np.loadtxt('data/hatomdat.txt', dtype=np.str)

for i in range (0, len(wasserstoff)):
    f, A = np.genfromtxt(wasserstoff[i], unpack=True)
    plt.plot()
    #plt.plot(f, A, 'ro', mew=0.5)
    plt.plot(f, A, 'b-', linewidth=0.5, label='Gemessenes Frequenzspektrum')
    plt.xlabel(r'$f/$Hz')
    plt.ylabel(r'$A$')
    plt.xlim()
    plt.ylim()
    plt.tight_layout()
    plt.legend()
    plt.grid()
    plt.savefig('build/hatom'+str(i)+'.pdf')
    plt.clf()
