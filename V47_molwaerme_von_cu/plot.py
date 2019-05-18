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

theta_D = 345

params, covariance_matrix = optimize.curve_fit(debyeFunction, T, C_V)
errors = np.sqrt(np.diag(covariance_matrix))
theta_D = params[0]
sigma_theta_D = errors[0]
theta_D_ufloat = ufloat(theta_D, sigma_theta_D)
print("theta_D = ", theta_D_ufloat)

Tlin = np.linspace(1, 60, 50)
C_Vplot = debyeFunction(Tlin, theta_D)
plt.plot(T, C_V, 'rx', label='Messwerte')
plt.plot(Tlin, C_Vplot, 'b-', label='Ausgleichsfunktion')
plt.grid()
plt.xlabel(r'$T$/K')
plt.ylabel(r'$C_V/$(J/mol K)')
plt.savefig('build/test.pdf')
