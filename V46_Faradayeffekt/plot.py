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

#Einlesen Daten erste dotierte Probe
grad1_dotiert_1296, minuten1_dotiert_1296, grad2_dotiert_1296, minuten2_dotiert_1296, lamb = np.genfromtxt('data/dotiert_1296.txt', unpack=True)
L = 0.0051

#Berechnung der Winkel erste dotierte Probe
theta1_dotiert_1296 = grad1_dotiert_1296 + minuten1_dotiert_1296/60
theta2_dotiert_1296 = grad2_dotiert_1296 + minuten2_dotiert_1296/60
#Winkel in rad umrechnen
theta1_dotiert_1296 = theta1_dotiert_1296/360 * 2*np.pi
theta2_dotiert_1296 = theta2_dotiert_1296/360 * 2*np.pi
#Berechnung des Drehwinkels
theta_dotiert_1296 = 1/2*(theta2_dotiert_1296 - theta1_dotiert_1296)
print(360*theta_dotiert_1296/(2*np.pi))
