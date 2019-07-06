import numpy as np
from scipy import integrate
LCDMf = lambda x: 1.0/np.sqrt(0.3*(1+x)**3+0.7)
np.vectorize(LCDMf)

def LCDMfint(z):
    return integrate.quad(LCDMf, 0, z)

LCDMfint=np.vectorize(LCDMfint)
z=np.arange(0,100)

an=LCDMfint(z)
print(an[0])