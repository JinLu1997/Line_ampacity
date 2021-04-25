### This file defines a function to calculate the line rating given conductor temperature Tc and other parameters

from math import sqrt, cos, acos, sin
# Returns the current required for a conductor temperature, Tc, and other parameters
def get_i(Tc, R_T_high, R_T_low, T_high, T_low, Ta, rho_f, D, epsilon, alpha, Q_se, theta, area):

    R_Tc = ((R_T_high - R_T_low) / (T_high - T_low)) * (Tc - T_low) + R_T_low

    # Conductor heat loss from natural convection, assuming no wind
    q_cn = 0.0205 * rho_f**0.5 * D**0.75 * (Tc - Ta)**1.25

    # Conductor heat loss from radiation
    q_r = 0.0178 * D * epsilon * (((Tc + 273)/100)**4 - ((Ta + 273)/100)**4)

    # Solar heat gain in conductor
    q_s = alpha * Q_se * sin(theta) * area

    I = sqrt((q_cn + q_r - q_s) / R_Tc)

    return I

### test data
D = 22.8                # Conductor diameter (mm)
area = D/1000           # Projected area of conducter per unit length (m^2/m)
rho_f = 1.029           # Density of air (kg/m^3)
H_e = 25                # Elevation of conductor above sea level (m)
Ta = 25.0               # Ambient temperature (C)
epsilon = 0.5           # Emissivity
alpha = 0.5             # Solar absorptivity
H_c = 72.5              # Altitude of sun (degrees)
Z_c = 139               # Azimuth of sun (degrees)
Z_l = 90.0              # Azimuth of electrical line (degrees): 90 degrees (or 270 degrees) for east-west
R_T_high = 8.688e-5     # Conductor unit resistance at high temperature reference (ohm/m)
R_T_low = 7.283e-5      # Conductor unit resistance at low temperature reference (ohm/m)
T_high = 75.0           # High temperature reference (C)
T_low = 25.0            # Low temperature reference (C)
Q_s = -42.2391 + 63.8044*H_c - 1.9220*H_c**2 + 3.46921e-2*H_c**3 - 3.61118e-4*H_c**4 + 1.94318e-6*H_c**5 - 4.07608e-9*H_c**6
K_solar = 1 + 1.148e-4*H_e - 1.108e-8*H_e**2
Q_se = K_solar * Q_s
theta = acos(cos(H_c) * cos(Z_c - Z_l))

### run get_i function with test data
# Given the max line temperature Tc (c)
Tc = 90
line_ampacity = get_i(Tc, R_T_high, R_T_low, T_high, T_low, Ta, rho_f, D, epsilon, alpha, Q_se, theta, area)
print('The line ampacity is: %.2f (A)' % line_ampacity)