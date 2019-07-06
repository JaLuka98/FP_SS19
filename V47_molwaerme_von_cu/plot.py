import numpy as np
from scipy import constants
import scipy.integrate as integrate
import uncertainties as unc
import scipy
import uncertainties.unumpy as unp
from scipy import optimize
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy import stats
from uncertainties import correlated_values
from matrix2latex import matrix2latex


def debyeFunction(T, theta_D):
    return [9*scipy.constants.R * (Temp/theta_D)**3 * integrate.quad(lambda x: x**4*np.exp(x)/(np.exp(x)-1)**2, 0.0001, theta_D/Temp)[0] for Temp in T]
    # [0] to only the return the nominal value, [1] would be the uncertainty


T = np.array([10, 16, 20, 25, 30, 35, 40, 50])
C_V = np.array([0.0555, 0.225, 0.462, 0.963, 1.693, 2.64, 3.74, 6.15])

#theta_D = 345
#Temperaturen und so berechnen

n=5.382 #mol. Folgt aus n=m/M mit der Masse der Probe und der molaren Masse von Cu

Rp, Rz=np.genfromtxt('data/data.txt', unpack=True)
Rp = unp.uarray(Rp, 0.1)
Rz = unp.uarray(Rz, 0.1)
I, t, U=np.genfromtxt('data/data2.txt', unpack=True)
I=unp.uarray(I, 0.1)
t=unp.uarray(t, 3)
U=unp.uarray(U, 0.1)
I41=np.append(0, I)
t41=np.append(0, t)
U41=np.append(0, U)
hr = ['$R_p/\Omega$', '$R_z/\Omega$','$I$/µA', '$t$/s', '$U/$V']
m = np.zeros((41, 10))
m[:,0] = unp.nominal_values(Rp)
m[:,1] = unp.std_devs(Rp)
m[:,2] = unp.nominal_values(Rz)
m[:,3] = unp.std_devs(Rz)
m[:,4] = unp.nominal_values(I41)
m[:,5] = unp.std_devs(I41)
m[:,6] = unp.nominal_values(t41)
m[:,7] = unp.std_devs(t41)
m[:,8] = unp.nominal_values(U41)
m[:,9] = unp.std_devs(U41)
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

Tp=0.00134*Rp**2+2.296*Rp -243.02 +273.15
Tz=0.00134*Rz**2+2.296*Rz -243.02 +273.15
I, t, U=np.genfromtxt('data/data2.txt', unpack=True)
I=unp.uarray(I, 0.1)
I*=1e-3 #jetzt in Ampere
t=unp.uarray(t, 3)
U=unp.uarray(U, 0.1)

E=I*U*t

DeltaT=[]
for i in range(1, 41):
    delta=Tp[i]-Tp[i-1]
    DeltaT.append(delta)


C_p=E/(n*np.asarray(DeltaT))

C_p41=np.append(0,C_p)
DeltaT41=np.append(0,DeltaT)
E41=np.append(0,E)

hr = ['$T_p/$K', '$\Delta T_p/$K','$E$/J', '$C_p$/\frac{mol}{kg K}']
m = np.zeros((41, 8))
m[:,0] = unp.nominal_values(Tp)
m[:,1] = unp.std_devs(Tp)
m[:,2] = unp.nominal_values(DeltaT41)
m[:,3] = unp.std_devs(DeltaT41)
m[:,4] = unp.nominal_values(E41)
m[:,5] = unp.std_devs(E41)
m[:,6] = unp.nominal_values(C_p41)
m[:,7] = unp.std_devs(C_p41)
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)


#C_v berechnen
def alphafunktion(T, a, b):
    return a/T+b

T, alpha=np.genfromtxt('data/alpha.txt', unpack=True)
alpha*=1e-6

V0=7.11*1e-6 #m3/mol
kappa=137.8*1e9 #N/m^2

params, covariance_matrix = optimize.curve_fit(alphafunktion, T, alpha)
a, b = correlated_values(params, covariance_matrix)
print('Fit für alpha:')
print('a=', a)
print('b=', b)

linspace=np.linspace(50, 320, 1000)
plt.plot(linspace, 1e6*alphafunktion(linspace, *params), 'b-', label='Ausgleichrechnung', linewidth=0.5)
plt.plot(T, alpha*1e6, 'rx', mew=0.5, label='Gegebene Werte')
plt.xlabel(r'$T/$K')
plt.ylabel(r'$\alpha/10^{-6}$K')
plt.axis([50,320,0,20])
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/alpha.pdf')
plt.clf()

alpha_T=a/Tp+b

C_v=C_p41-alpha_T**2*kappa*V0*Tp

hr = ['$T$/K', '$\alpha$/10^-6 °']
m = np.zeros((24, 2))
m[:,0] = T
m[:,1] = alpha*1e6
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

hr = ['$C_V$/J/kg K']
m = np.zeros((41, 2))
m[:,0] = unp.nominal_values(C_v)
m[:,1] = unp.std_devs(C_v)
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)


plt.errorbar(unp.nominal_values(Tp[1:41]), unp.nominal_values(C_v[1:41]), yerr=unp.std_devs(C_v[1:41]), xerr=unp.std_devs(Tp[1:41]), fmt='rx',mew=0.5, ecolor='b', elinewidth=0.5,  label='Errechnete Werte') #Irgendwas klappt hier mit den errorbars noch nicht so richtig
plt.xlabel(r'$T/$K')
plt.ylabel(r'$C_v$ in J/(mol K)')
#plt.axis([50,320,0,0.00002])
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/cv.pdf')
plt.clf()



#Dieses komische Debye Zeug

params, covariance_matrix = optimize.curve_fit(debyeFunction, unp.nominal_values(Tp[1:21]), unp.nominal_values(C_v[1:21]), sigma=unp.std_devs(C_v[1:21]), absolute_sigma=True)
errors = np.sqrt(np.diag(covariance_matrix))
theta_D = params[0]
sigma_theta_D = errors[0]
theta_D_ufloat = ufloat(theta_D, sigma_theta_D)
print("theta_D = ", theta_D_ufloat)

Tlin = np.linspace(-1, 180, 200)
C_Vplot = debyeFunction(Tlin, theta_D)
plt.errorbar(unp.nominal_values(Tp[1:21]), unp.nominal_values(C_v[1:21]), label='Messwerte', yerr=unp.std_devs(C_v[1:21]), fmt='rx',mew=0.5, ecolor='b', elinewidth=0.5)
plt.plot(Tlin, C_Vplot, 'g-', label='Ausgleichsfunktion', linewidth=0.5)
plt.grid()
plt.legend()
plt.axis([0, 180, -1, 25])
plt.xlabel(r'$T$/K')
plt.ylabel(r'$C_V/$(J/mol K)')
plt.savefig('build/debye.pdf')
plt.clf()

# Calculating Debye Frequency and Temperature analytically
N_A = 6.02214179*1e23 # in 1/mol
M = 63.55*1e-3 # in kg/mol
rho = 8960 # in kg/m^3
number_density = N_A*rho/M

v_long = 4.7*1e3 # in m/s
v_trans = 2.26*1e3 # in m/s

omega_D = ((18*np.pi**2*number_density) * (1/v_long**3+2/v_trans**3)**(-1))**(1/3)
print('Debye-Frequenz', omega_D)

k_B = 1.3806504*1e-23
hbar = 1.054571628*1e-34

Theta_D = hbar*omega_D/k_B
print('Debye-Temperatur', Theta_D)
