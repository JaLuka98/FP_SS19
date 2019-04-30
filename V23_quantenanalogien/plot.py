import numpy as np
import scipy.optimize
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy import optimize
import matplotlib.pyplot as plt
from uncertainties import ufloat
from uncertainties import correlated_values
from matrix2latex import matrix2latex


#Polarplot

alpha, A1, A2, A3=np.genfromtxt('data/peaks.txt', unpack=True)

hr = ['$\alpha$/°', '$A_1$', '$A_2$','$A_3$']
m = np.zeros((19, 4))
m[:,0] = alpha
m[:,1] = A1
m[:,2] = A2
m[:,3] = A3
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

#Vermutlich P_1,0
plt.figure(3)
phirange = np.linspace(0, 2*np.pi, 1000)
plt.polar(np.arccos(1/2*np.cos(alpha/360*2*np.pi)-1/2), A1/np.max(A1), "rx", mew=0.5, label="Messwerte")
plt.polar(phirange, np.abs(np.cos(phirange))/max(np.abs(np.cos(phirange))), "b-", label="Theorie", linewidth=0.5)
plt.legend()
plt.savefig("build/polar1.pdf")
plt.clf()

#Vermutlich P_2,0
plt.figure(3)
plt.polar(np.arccos(1/2*np.cos(alpha/360*2*np.pi)-1/2), A2/np.max(A2), "rx", mew=0.5, label="Messwerte")
plt.polar(phirange, 1/2*(3*(np.abs(np.cos(phirange)))**2-1)/max(1/2*(3*(np.abs(np.cos(phirange)))**2-1)), "b-", label="Theorie", linewidth=0.5)
plt.legend()
plt.savefig("build/polar2.pdf")
plt.clf()

#Vermutlich P_3,0... der passt mit dem gemogelten Minus und ein bisschen Phantasie
plt.figure(3)
plt.polar(np.arccos(1/2*np.cos(alpha/360*2*np.pi)-1/2), A3/np.max(A3), "rx", mew=0.5, label="Messwerte")
plt.polar(phirange, (np.abs(np.cos(phirange))/2*(5*np.cos(phirange)**2-3))/max((np.cos(phirange))/2*(5*((np.cos(phirange)))**2-3)), "b-", label="Theorie", linewidth=0.5)
#Das Minus gehört da eigentlich nicht hin, aber sonst kommt da nichts raus was auch
#nur annähernd zu den Messwerten passt...
plt.legend()
plt.savefig("build/polar3.pdf")
plt.clf()


# Zylinderkette
def linfit(x, a, b):
    return a*x+b


zylinder, f1, a1, b1, A1, f2, a2, b2, A2 = np.genfromtxt('data/roehre.txt', unpack=True)

deltaf = 1000*(f2-f1)
zylinderanzahl = np.log(zylinder)
L = np.linspace(5,60,12)
einsDurchZweiL = 1/(2*L)
Llin = np.linspace(5,60,1000)
#linspace = np.linspace(-0.1, 2.7, 500)
params, covariance_matrix = optimize.curve_fit(linfit, einsDurchZweiL, deltaf)
a, b = correlated_values(params, covariance_matrix)
print('Fit zur Schallgeschwindigkeitsbestimmung')
print('a=', a)
print('b=', b)

plt.plot(einsDurchZweiL, deltaf/1000, 'rx', mew=0.5, label='Messwerte') # Durch 1000 für Kilohertz
plt.plot(1/(2*Llin), linfit(1/(2*Llin), *params)/1000, 'b-', label='Ausgleichrechnung', linewidth=0.5)
#plt.yscale('log')
#plt.xscale('log')
plt.xlabel(r'$\frac{1}{2L}\cdot$m')
plt.ylabel(r'$\Delta f$/kHz')
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/roehre.pdf')
plt.clf()

hr = ['Zylinderanzahl', 'Resonatorlänge/cm' '$f_1/$Hz', '$f_2/$Hz', '$\Delta f/$Hz']
m = np.zeros((12, 5))
m[:,0] = zylinder
m[:,1] = 5*zylinder
m[:,2] = f1*1000
m[:,3] = f2*1000
m[:,4] = 1000*(f2-f1)
t = matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

zylinder = np.genfromtxt('data/zyldat.txt', dtype=np.str)

for i in range(0, len(zylinder)):
    f, A = np.genfromtxt(zylinder[i], unpack=True)
    #plt.plot(f, A, 'ro', mew=0.5)
    plt.plot(f, A, 'b-', linewidth=0.5, label='Frequenzspektrum')
    f_osz = np.genfromtxt('data/oszispektrum.txt', usecols=(i))
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

f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12 = np.genfromtxt('data/oszispektrum.txt', unpack=True)
hr = ['$f_1$/Hz','$f_2$/Hz','$f_3$/Hz','$f_4$/Hz','$f_5$/Hz','$f_6$/Hz','$f_7$/Hz',
'$f_8$/Hz','$f_9$/Hz','$f_10$/Hz','$f_11$/Hz','$f_12$/Hz']
m = np.zeros((12, 13))
m[:,0] = f1*1000
m[:,1] = f2*1000
m[:,2] = f3*1000
m[:,3] = f4*1000
m[:,4] = f5*1000
m[:,5] = f6*1000
m[:,6] = f7*1000
m[:,7] = f8*1000
m[:,8] = f9*1000
m[:,9] = f10*1000
m[:,10] = f11*1000
m[:,11] = f12*1000
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)


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

hr = ['$f_\symup{res}$/Hz']
m = np.zeros((10, 1))
m[:,0] = f_osz*1000
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)


wasserstoff=np.genfromtxt('data/hatomdat.txt', dtype=np.str)

for i in range (0, len(wasserstoff)):
    f, A = np.genfromtxt(wasserstoff[i], unpack=True)
    plt.plot()
    #plt.plot(f, A, 'ro', mew=0.5)
    plt.plot(f, A, 'b-', linewidth=0.5, label='Frequenzspektrum')
    plt.xlabel(r'$f/$Hz')
    plt.ylabel(r'$A$')
    plt.xlim()
    plt.ylim()
    plt.tight_layout()
    plt.legend()
    plt.grid()
    plt.savefig('build/hatom'+str(i)+'.pdf')
    plt.clf()
