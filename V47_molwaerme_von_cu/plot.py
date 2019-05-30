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
I, t, U=np.genfromtxt('data/data2.txt', unpack=True)
I41=np.append(I, 0)
t41=np.append(t, 0)
U41=np.append(U, 0)
hr = ['$R_p/\Omega$', '$R_z/\Omega$','$I$/µA', '$t$/s', '$U/$V']
m = np.zeros((41, 5))
m[:,0] = Rp
m[:,1] = Rz
m[:,2] = I41
m[:,3] = t41
m[:,4] = U41
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

Tp=0.00134*Rp**2+2.296*Rp -243.02 +273.15
Tz=0.00134*Rz**2+2.296*Rz -243.02 +273.15
I, t, U=np.genfromtxt('data/data2.txt', unpack=True)
I*=1e-3 #jetzt in Ampere
E=I*U*t

DeltaT=[]
for i in range(1, 41):
    delta=Tp[i]-Tp[i-1]
    DeltaT.append(delta)


C_p=E/(n*np.asarray(DeltaT))

C_p41=np.append(C_p, 0)
DeltaT41=np.append(DeltaT, 0)
E41=np.append(E, 0)

hr = ['$T_p/$K', '$\Delta T_p/$K','$E$/J', '$C_p$/\frac{mol}{kg K}']
m = np.zeros((41, 4))
m[:,0] = Tp
m[:,1] = DeltaT41
m[:,2] = E41
m[:,3] = C_p41
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

hr = ['$T$/K', '$\alpha$/10^-6 °'] #hier vielleicht noch die werte für C_v rein?
m = np.zeros((24, 2))
m[:,0] = T
m[:,1] = alpha*1e6
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

hr = ['$C_V$/J/kg K']
m = np.zeros((41, 2))
m[:,0] = unp.nominal_values(C_v)
m[:,1] = unp.std_devs(C_v)
t=matrix2latex(m, headerRow=hr, format='%.5f')
print(t)

plt.errorbar(Tp, unp.nominal_values(C_v), fmt='rx',mew=0.5, label='Errechnete Werte') #Irgendwas klappt hier mit den errorbars noch nicht so richtig
plt.xlabel(r'$T/$K')
plt.ylabel(r'$C_v$ in J/(mol K)')
#plt.axis([50,320,0,0.00002])
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/cv.pdf')
plt.clf()



#Dieses komische Debye Zeug

params, covariance_matrix = optimize.curve_fit(debyeFunction, Tp[0:20], unp.nominal_values(C_v[0:20]))
errors = np.sqrt(np.diag(covariance_matrix))
theta_D = params[0]
sigma_theta_D = errors[0]
theta_D_ufloat = ufloat(theta_D, sigma_theta_D)
print("theta_D = ", theta_D_ufloat)

Tlin = np.linspace(-1, 180, 200)
C_Vplot = debyeFunction(Tlin, theta_D)
plt.plot(Tp[0:20], unp.nominal_values(C_v[0:20]), 'rx', mew=0.5, label='Messwerte')
plt.plot(Tlin, C_Vplot, 'b-', label='Ausgleichsfunktion', linewidth=0.5)
plt.grid()
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
