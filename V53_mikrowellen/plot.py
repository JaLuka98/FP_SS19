import numpy as np
import scipy.optimize
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy import optimize
import matplotlib.pyplot as plt
from uncertainties import ufloat


def parabola(x, a, b, c):
    return a*x*x + b*x + c


#Frequenzmessung
f_dip=9006 #MHz

#3dB Mehtode
d1=69.6 #mm
d2=67.9

#Abschwächer Methode
A_3=42  #dB

f, U_l, U_max, U_r, A = np.genfromtxt('data/moden.txt', unpack=True, comments='#')

jet= plt.get_cmap('jet')
colors = iter(jet(np.linspace(0,1,10)))

for i in range(0, 3):
    x = np.array([U_r[i], U_max[i], U_l[i]])
    xlin = np.linspace(U_r[i], U_l[i], 1000)
    y = np.array([0, A[i], 0])
    print(x)
    print(y)
    params, covariance_matrix = optimize.curve_fit(parabola, x, y)
    # errors = np.sqrt(np.diag(covariance_matrix)) because 3 points fit perfectly
    print("a, b, c = %.1f" % params[i])
    plt.plot(x, y, 'x', color=next(colors))
    plt.plot(xlin, parabola(xlin, *params), color=next(colors))
    plt.grid()
    plt.xlabel(r'$U_\mathrm{refl}$/V')
    plt.ylabel(r'$A/$mV')
    plt.savefig('build/modes.pdf')
