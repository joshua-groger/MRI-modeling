# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 16:50:14 2020

@author: joshu
"""

import numpy as np
from numpy import pi, log10, sqrt
import matplotlib.pyplot as plt

#

# Physical Constants
mu_0 = 4 * pi * 10 ** -7  # H/m
copResis = 1.69e-8  # Ohm*m
bw = 10e6  # Hertz

# Input Coil Parameters

coilDia = 0.113                                                             # m
wireDia = 4.02e-4                                                              # m
omega   = 2*pi*2.98e8                                                               # rad/s - resonant radial frequency
netImp  = 50                                                                # Ohm - impedance to match

width   = 0.03                                                              # m

#%%

# Derived values

# These values are calculated using formulas in Barral's paper
# There are equations available for both circular and rectangular coils

# Coil properties
copSkin = sqrt(2 * copResis / (mu_0 * omega))  # m - AC skin depth

# for wire
induc = mu_0 * coilDia / 2 * (log10(8 * coilDia / wireDia) - 2)             # H - coil inductance
# resis = coilDia * copResis / (wireDia * copSkin)                            # Ohm - coil resistance

# # for foil
resis = pi * coilDia * copResis / (2 * width * copSkin)

# Resonant capacitance
#
# The basic formula for a resonant LC circuit is omega_resonant = 1 / sqrt(inductance * capacitance)
# The capacitance used here is an offshoot of this to put us in the right direction to network matching.
# If we take C_res and add (or remove) a capacitor in series with the coil, we can move the impedance to
# where it lies on the unit admitance circle on a smith chart. This is where Re{1/Y} = 1.
# Adding another capacitor in parallel with the coil removes the reactance component from the impedance,
# giving us a "perfect" network match.
# By some miracle this still resonates at the desired frequency, but with the desired imput impedance.
#
# 1.) C_res = 1/(omega_res^2 * inductance)
# 2.) "Remove" capacitance to hit the circle - Z = 1/jwC_res + R + jwL + 1/jwC_ser
#       Y = 1/Z = 1 / (1/jwC_res + R + jwL + 1/jwC_ser)
#       Choose C_ser such that Re{Y} = Yo = 1/Zo
#       Re{Y} = Re{1 / (R + jX)} = R / (R^2 + X^2)
#       With some algebra R*Zo = R^2 + X^2
#       Solve for X = sqrt(R*Zo - R^2) = -1/wC_res + wL - 1/wC_ser
#       Solve for C_ser = 1 / (w * (wL - 1/wC_res - sqrt(R*Zo - R^2)))
# 3.) Replace C_res and C_ser with C_eq = (C_ser * C_res) / (C_ser + C_res)
#       This simplifies to C_eq = 1 / (w * (w * L - sqrt(R * Zo - R**2)))
# 4.) Add a parallel capacitor to remove reactance
#       0 = Im{Y + 1/jwC_par}
#       0 = Im{1/(R+jX) + jB} = -X / (R^2 + X^2) + B
#       B = X / (R^2 + X^2) = wC_par
#       C_par = (wL - 1/wC_eq) / (w * (R^2 + (wL - 1/wC_eq)^2))

# capRes = 1 / (omega**2 * induc)
capEq = 1 / (omega * (omega * induc - sqrt(resis * netImp - resis ** 2)))  # F

# Parallel capacitance
capPar = (omega * induc - 1 / (omega * capEq)) / (omega * (resis ** 2 + (omega * induc - 1 / (omega * capEq)) ** 2))


# Input impedance
# This doesn't need to be a function necessarily, but it computes the input impedance for the loop
def impedance(omega, resistance, inductance, capacitance, trim):
    return 1 / (1 / (resistance + 1j * omega * inductance + 1 / (1j * omega * capacitance)) + 1j * omega * trim)


#

# Display Impedance Spectrum

# I calculate the impedance for a range of frequencies around our desired frequency
omegas = np.linspace(omega - bw / 2, omega + bw / 2)
impeds = np.asarray([impedance(w, resis, induc, capEq, capPar) for w in omegas])

# Then I plot it
plt.plot(omegas, impeds.real, label='Resistance')
plt.plot(omegas, impeds.imag, '--', label='Reactance')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('Resistance (Ohms)')
plt.legend()

print(f'Resistance: {resis}\nInductance: {induc}\nResonant Capacitor: {capEq}\nParallel Capacitor: {capPar}')
