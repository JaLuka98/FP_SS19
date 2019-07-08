﻿import numpy as np
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


def linearFunction(x,a):
    return a*x


#Einlesen Daten reine Probe
grad1_rein, minuten1_rein, grad2_rein, minuten2_rein, lamb = np.genfromtxt('data/rein.txt', unpack=True)
L = 0.0051

#Berechnung der Winkel reine Probe
theta1_rein = grad1_rein + minuten1_rein/60
theta2_rein = grad2_rein + minuten2_rein/60
#Winkel in rad umrechnen
theta1_rein = theta1_rein/360 * 2*np.pi
theta2_rein = theta2_rein/360 * 2*np.pi
#Berechnung des Drehwinkels
theta_rein = 1/2*(theta2_rein - theta1_rein)
print(360*theta_rein/(2*np.pi))

hr = [r'$\lambda/\mu$m', r'$\theta_1$/rad', r'$\theta_2$/rad', r'$\theta$/rad']
size = np.size(lamb)
print(lamb)
m = np.zeros((size, 4))
m[:, 0] = lamb
m[:, 1] = theta1_rein[0:size]
m[:, 2] = theta2_rein[0:size]
m[:, 3] = theta_rein[0:size]
t = matrix2latex(m, headerRow=hr, format='%.4f')
print(t)

theta_frei_rein = theta_rein/L  #  in rad/m
plt.plot(lamb, theta_frei_rein, 'rx')
plt.xlabel(r'$\lambda/\mu$m')
plt.ylabel(r'$\frac{\theta}{L}/\frac{\mathrm{rad}}{\mathrm{mm}}$')
plt.grid()
plt.axis([0.875, 2.875, -2.5, 32.5])
plt.savefig('build/rein.pdf')
plt.clf()

#Einlesen Daten erste dotierte Probe
grad1_dotiert_136, minuten1_dotiert_136, grad2_dotiert_136, minuten2_dotiert_136, lamb = np.genfromtxt('data/dotiert_136.txt', unpack=True)
L = 0.00136

#Berechnung der Winkel erste dotierte Probe
theta1_dotiert_136 = grad1_dotiert_136 + minuten1_dotiert_136/60
theta2_dotiert_136 = grad2_dotiert_136 + minuten2_dotiert_136/60
#Winkel in rad umrechnen
theta1_dotiert_136 = theta1_dotiert_136/360 * 2*np.pi
theta2_dotiert_136 = theta2_dotiert_136/360 * 2*np.pi
#Berechnung des Drehwinkels
theta_dotiert_136 = 1/2*(theta2_dotiert_136 - theta1_dotiert_136)
print(360*theta_dotiert_136/(2*np.pi))

hr = [r'$\lambda/\mu$m', r'$\theta_1$/rad', r'$\theta_2$/rad', r'$\theta$/rad']
size = np.size(lamb)
print(lamb)
m = np.zeros((size, 4))
m[:, 0] = lamb
m[:, 1] = theta1_dotiert_136[0:size]
m[:, 2] = theta2_dotiert_136[0:size]
m[:, 3] = theta_dotiert_136[0:size]
t = matrix2latex(m, headerRow=hr, format='%.4f')
print(t)

theta_frei_dotiert_136 = theta_dotiert_136/L  #  in rad/m
plt.plot(lamb, theta_frei_dotiert_136, 'rx')
plt.xlabel(r'$\lambda/\mu$m')
plt.ylabel(r'$\frac{\theta}{L}/\frac{\mathrm{rad}}{\mathrm{mm}}$')
plt.grid()
plt.axis([0.875, 2.875, 15, 50])
plt.savefig('build/dotiert_136.pdf')
plt.clf()

#Einlesen Daten zweite dotierte Probe
grad1_dotiert_1296, minuten1_dotiert_1296, grad2_dotiert_1296, minuten2_dotiert_1296, lamb = np.genfromtxt('data/dotiert_1296.txt', unpack=True)
L = 0.001296

#Berechnung der Winkel zweite dotierte Probe
theta1_dotiert_1296 = grad1_dotiert_1296 + minuten1_dotiert_1296/60
theta2_dotiert_1296 = grad2_dotiert_1296 + minuten2_dotiert_1296/60
#Winkel in rad umrechnen
theta1_dotiert_1296 = theta1_dotiert_1296/360 * 2*np.pi
theta2_dotiert_1296 = theta2_dotiert_1296/360 * 2*np.pi
#Berechnung des Drehwinkels
theta_dotiert_1296 = 1/2*(theta2_dotiert_1296 - theta1_dotiert_1296)
print(360*theta_dotiert_1296/(2*np.pi))

hr = [r'$\lambda/\mu$m', r'$\theta_1$/rad', r'$\theta_2$/rad', r'$\theta$/rad']
size = np.size(lamb)
print(lamb)
m = np.zeros((size, 4))
m[:, 0] = lamb
m[:, 1] = theta1_dotiert_1296[0:size]
m[:, 2] = theta2_dotiert_1296[0:size]
m[:, 3] = theta_dotiert_1296[0:size]
t = matrix2latex(m, headerRow=hr, format='%.4f')
print(t)

theta_frei_dotiert_1296 = theta_dotiert_1296/L  #  in rad/m
plt.plot(lamb, theta_frei_dotiert_1296, 'rx')
plt.xlabel(r'$\lambda/\mu$m')
plt.ylabel(r'$\frac{\theta}{L}/\frac{\mathrm{rad}}{\mathrm{mm}}$')
plt.grid()
plt.axis([0.875, 2.875, 25, 95])
plt.savefig('build/dotiert_1296.pdf')
plt.clf()

#######################################
### Bestimmung der effektiven Masse ###
#######################################

theta_136 = theta_frei_dotiert_136 - theta_frei_rein
params, covariance_matrix = optimize.curve_fit(linearFunction, lamb**2, theta_136)
a = correlated_values(params, covariance_matrix)
print('------------------------------------------------------------')
print('Fit parameter für die erste Probe (136): ', a)
print('------------------------------------------------------------')

linLin = np.linspace(0, 10, 1000)
plt.plot(lamb**2, theta_136, 'rx', label='Messwerte', zorder=2)
plt.plot(linLin, linearFunction(linLin, *params), 'b-', label='Anpassungsfunktion', zorder=3)
plt.legend()
plt.grid()
plt.axis([0, 8, 0, 50])
#plt.gcf().subplots_adjust(bottom=0.18)
plt.xlabel(r'$\lambda^2 / \mu m^2$')
plt.ylabel(r'$\frac{\Theta}{L} / \frac{\mathrm{rad}}{\mathrm{m}}$')
plt.savefig('build/differenz_136.pdf')
plt.clf()

theta_1296 = theta_frei_dotiert_1296 - theta_frei_rein
params, covariance_matrix = optimize.curve_fit(linearFunction, lamb**2, theta_1296)
a = correlated_values(params, covariance_matrix)
print('------------------------------------------------------------')
print('Fit parameter für die erste Probe (136): ', a)
print('------------------------------------------------------------')

linLin = np.linspace(0, 10, 1000)
plt.plot(lamb**2, theta_1296, 'rx', label='Messwerte', zorder=2)
plt.plot(linLin, linearFunction(linLin, *params), 'b-', label='Anpassungsfunktion', zorder=3)
plt.legend()
plt.grid()
plt.axis([0, 8, 0, 100])
#plt.gcf().subplots_adjust(bottom=0.18)
plt.xlabel(r'$\lambda^2 / \mu m^2$')
plt.ylabel(r'$\frac{\Theta}{L} / \frac{\mathrm{rad}}{\mathrm{m}}$')
plt.savefig('build/differenz_1296.pdf')
plt.clf()
