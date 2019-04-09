import numpy as np
import scipy.optimize
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy import optimize
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy import stats


def parabola(x, a, b, c):
    return a*x*x + b*x + c


#Frequenzmessung
f_dip=9006 #MHz

#3dB Mehtode
d1=69.6 #mm
d2=67.9

#Abschwächer Methode
A_3=42  #dB

# Untersuchung der Moden

f, U_l, U_max, U_r, A = np.genfromtxt('data/moden.txt', unpack=True, comments='#')

jet = plt.get_cmap('jet')
colors = iter(jet(np.linspace(0, 1, 10)))

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
    plt.clf()

# Frequenz- und Wellenlängen untersuchen

# Abmessungen des Hohlleiters hier einfügen (in millimetern)
a = 22.86
b = 10

x = np.genfromtxt('data/wellenlaenge.txt', unpack=True, comments='#')
x1 = x[1] - x[0]
x2 = x[2] - x[1]
x3 = x[3] - x[2]
x = np.array([x1, x2, x3])  # x in mm
lam_g = ufloat(np.mean(x), stats.sem(x))  # default ddof = 1
print('Wellenlänge im Hohlleiter', lam_g)
lam_c = 2*a
print('Grenzwellenlänge des Hohlleiters', lam_c)
lam_0 = 1/unp.sqrt(1/lam_g**2 + 1/(2*a)**2)
print('Wellenlänge im freien Raum', lam_0)
f_exp = 3*1e11*unp.sqrt(1/lam_g**2 + 1/(2*a)**2)
print('f_exp ist', f_exp)
v_ph = f_exp*lam_g
print('v_phase', v_ph)

# Betrachtung der Dämpfung

dB_selbst, d_selbst = np.genfromtxt('data/daempfung.txt', unpack=True, comments='#')
d_dort, dB_dort = np.genfromtxt('data/daempfungbild.txt', unpack=True, comments='#')

plt.plot(d_selbst, dB_selbst, 'rx', label='Eigene Messung')
plt.plot(d_dort, dB_dort, 'bx', label='Herstellerangabe')
plt.grid()
plt.legend()
plt.xlabel(r'$d/$mm')
plt.ylabel(r'$$dB')
plt.savefig('build/daempfung.pdf')
