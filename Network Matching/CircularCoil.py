# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 16:50:14 2020

@author: joshu
"""

import numpy as np
from numpy import pi, log10, sqrt
import matplotlib.pyplot as plt

#%%

# Constants
mu_0 = 4 * pi * 10**-7          # H/m
copResis = 1.69e-8              # Ohm*m
bw = 10e6                       # Hertz

# Parameters
coilDia = 0.1                   # m
wireDia = 1e-3                  # m
omega   = 1e8                   # rad/s
netImp  = 50                    # Ohm

# Coil properties
copSkin = sqrt(2 * copResis / (mu_0 * omega))

# for wire
induc = mu_0 * coilDia / 2 * (log10(8 * coilDia / wireDia) - 2)
resis = coilDia * copResis / (wireDia * copSkin)

# # for foil
# resis = pi * coilDia * copResis / (2 * width * copSkin)

# Resonant capacitance
capRes = 1 / (omega * (omega * induc - sqrt(resis * netImp - resis**2)))
#capRes = 1 / (omega**2 * induc)

# # Series capacitance
# capSer = sqrt(1 / (netImp * resis - resis**2)) / omega

# Parallel capacitance
capPar =  (omega * induc - 1 / (omega * capRes)) / (omega * (resis**2 + (omega * induc - 1 / (omega * capRes))**2))
#capPar = 10e-6

# Input impedance
def impedance(omega, resistance, inductance, capacitance, trim):
    return 1 / (1 / (resistance + 1j * omega * inductance + 1 / (1j * omega * capacitance)) + 1j * omega * trim)

omegas = np.linspace(omega-bw/2, omega+bw/2)
impeds = np.asarray([impedance(w, resis, induc, capRes, capPar) for w in omegas])

plt.plot(omegas, impeds.real, label='Resistance')
plt.plot(omegas, impeds.imag, '--', label='Reactance')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('Resistance (Ohms)')
plt.legend()