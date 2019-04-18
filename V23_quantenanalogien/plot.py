import numpy as np
import scipy.optimize
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy import optimize
import matplotlib.pyplot as plt
from uncertainties import ufloat


zylinder=np.loadtxt('data/zyldat.txt', dtype=np.str)

for i in range (0, len(zylinder)):
    f, A = np.genfromtxt(zylinder[i], unpack=True)
    plt.plot()
    #plt.plot(f, A, 'ro', mew=0.5)
    plt.plot(f, A, 'b-', linewidth=0.5)
    plt.xlabel(r'$f/$Hz')
    plt.ylabel(r'$A$')
    plt.xlim()
    plt.ylim()
    plt.tight_layout()
    plt.legend()
    plt.grid()
    plt.savefig('build/zyl'+str(i+1)+'.pdf')
    plt.clf()

wasserstoff=np.loadtxt('data/hatomdat.txt', dtype=np.str)

for i in range (0, len(wasserstoff)):
    f, A = np.genfromtxt(wasserstoff[i], unpack=True)
    plt.plot()
    #plt.plot(f, A, 'ro', mew=0.5)
    plt.plot(f, A, 'b-', linewidth=0.5)
    plt.xlabel(r'$f/$Hz')
    plt.ylabel(r'$A$')
    plt.xlim()
    plt.ylim()
    plt.tight_layout()
    plt.legend()
    plt.grid()
    plt.savefig('build/hatom'+str(i)+'.pdf')
    plt.clf()
