#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 08:40:12 2019

@author: christian
"""

import numpy as np
import scipy.optimize as optimize


def HW3_15(rho, V):
    Sw = 5500
    b = 196
    W = 636600
    CDo = 0.020
    CDoL = 0.0
    e = 0.72
    RA = b*b/Sw
    W_Sw = W/Sw
    rho_o = 0.0023769
    TA = 250000*(rho/rho_o)**0.8
    TA_W = TA/W
    # Eq. (3.4.8)
    # PA is TA*V and PR is given by Eq. (3.3.7)
    V_c = (TA_W*V - (CDo*rho*V*V*V)/(2.*(W_Sw)) - CDoL*V -
           (2.*W_Sw)/(np.pi*e*RA*rho*V))
    return V_c


def HW3_16(rho):
    Sw = 5500
    b = 196
    W = 636600
    CDo = 0.020
    CDoL = 0.0
    e = 0.72
    RA = b*b/Sw
    W_Sw = W/Sw
    rho_o = 0.0023769
    TA = 250000*(rho/rho_o)**0.8
    TA_W = TA/W
    # Found by taking the derivative of Eq. (3.4.8) with PA=TA*V and PR given
    # by Eq. (3.3.7) with respect to V

    def dVc_dV(V):
        der = (TA/W - (3.*CDo*rho*V*V)/(2.*W_Sw) - CDoL +
               (2.*W_Sw)/(np.pi*e*RA*rho*V*V))
        return abs(der)
    # Optimize to set this equal to zero and use the result in Eq. (3.4.8)
    V = optimize.minimize(dVc_dV, 200).x
    V_c = (TA_W*V - (CDo*rho*V*V*V)/(2.*(W_Sw)) - CDoL*V -
           (2.*W_Sw)/(np.pi*e*RA*rho*V))
    return V, V_c


def HW3_19(rho):
    RA = 12
    W = 200
    PA = 0.7*0.33*550
    PR = PA
    CDo = 0.010
    e = 0.80
    rho_o = 0.0023769
    # Solve Eq. (3.3.11) for Sw
    Sw = (W*W*W)/(PR*PR*rho_o)*(32.*CDo**0.5)/(3.*np.pi*e*RA)**1.5
    W_Sw = W/Sw
    # From Eq. (3.3.12)
    V_MDV = np.sqrt(W_Sw/rho_o)*(2.)/(12.*np.pi*e*RA*CDo)**0.25
    # Scale the answer to the given density
    P_input = 0.33*np.sqrt(rho_o/rho)
    return Sw, V_MDV, P_input


def HW3_20():
    W = 200
    # Assuming the same efficiency of mechanical work to power as 3.19
    PA = 0.7*0.6*550
    PR = 0.7*0.33*550
    V_c = (PA - PR)/W
    return V_c


def HW3_42(V_wind):
    Sw = 5500
    b = 196
    W = 636600
    W_Sw = W/Sw
    CDo = 0.02
    e = 0.72
    RA = b*b/Sw
    rho = 0.0023769
    R_Vhw = V_wind/np.sqrt(W_Sw/rho)

    def RVpoly(RV):
        res = (CDo*RV**5 - (3.*CDo*R_Vhw*RV**4)/2. -
               (4.*RV)/(np.pi*e*RA) + 2.*R_Vhw/(np.pi*e*RA))
        return abs(res)
    # Follow Example 3.7.1
    if V_wind == 0.:
        R_V = (4./(np.pi*e*RA*CDo))**.25
        R_Gmax = np.sqrt(np.pi*e*RA/(4.*CDo))
    if V_wind > 0.:
        R_V = optimize.minimize(RVpoly, 2.).x
        R_Gmax = ((1. - R_Vhw/R_V)*((CDo*R_V*R_V)/2. +
                                    2./(np.pi*e*RA*R_V*R_V))**-1)
    if V_wind < 0.:
        R_V = optimize.minimize(RVpoly, 2.).x
        R_Gmax = ((1. - R_Vhw/R_V)*((CDo*R_V*R_V)/2. +
                                    2./(np.pi*e*RA*R_V*R_V))**-1)
    V_BG = R_V*np.sqrt(W_Sw/rho)
    return V_BG, R_Gmax


def HW3_45(air_num):
    if air_num == 7:
        Sw = 950
        W = 73000
        CLmax = 1.2
        rho = 0.0023769
    if air_num == 10:
        Sw = 5500
        W = 636600
        CLmax = 1.2
        rho = 0.0023769
    if air_num == 17:
        Sw = 945
        W = 40000
        CLmax = 1.2
        rho = 0.0023769
    if air_num == 19:
        Sw = 775
        W = 200
        CLmax = 1.2
        rho = 0.0023769
    V_min = np.sqrt(2./CLmax)*np.sqrt(W/Sw/rho)
    print("3.%d: V_min=%.12f ft/s" % (air_num, V_min))


def HW3_47(rho):
    Sw = 5500
    b = 196
    W = 636600
    CDo = 0.020
    CDoL = 0.0
    e = 0.72
    RA = b*b/Sw
    W_Sw = W/Sw
    CLmax = 1.2
    npll = 4.0
    g = 32.2
    V_M = np.sqrt(2.*npll/CLmax)*np.sqrt(W_Sw/rho)
    R_min = 2.*npll/(CLmax*np.sqrt(npll*npll - 1.))*W_Sw/(rho*g)
    Om_max = np.sqrt((CLmax/2.)*(npll - 1./npll))*g*np.sqrt(rho/W_Sw)
    print("V_M: %.12f\nR_min: %.12f\nOmega_max: %.12f" % (V_M, R_min, Om_max))


print("HW 3.15")
print("------------------------")
print("Sea Level")
print("V=400 ft/s", "V_c=", HW3_15(0.0023769, 400), "ft/s")
print("V=600 ft/s", "V_c=", HW3_15(0.0023769, 600), "ft/s")
print("V=800 ft/s", "V_c=", HW3_15(0.0023769, 800), "ft/s")
print("30,000 ft")
print("V=400 ft/s->", "V_c=", HW3_15(0.00089068, 400), "ft/s")
print("V=600 ft/s->", "V_c=", HW3_15(0.00089068, 600), "ft/s")
print("V=800 ft/s->", "V_c=", HW3_15(0.00089068, 800), "ft/s")
print("------------------------\n")

print("HW 3.16")
print("------------------------")
print("Sea Level")
print("V=%.12f ft/s\nV_cmax=%.12f ft/s\t(%.12f fpm)" % (HW3_16(0.0023769)[0],
                                                        HW3_16(0.0023769)[1],
                                                        HW3_16(0.0023769)[1]*60))
print("30,000 ft")
print("V=%.12f ft/s\nV_cmax=%.12f ft/s\t(%.12f fpm)" % (HW3_16(.00089068)[0],
                                                        HW3_16(.00089068)[1],
                                                        HW3_16(.00089068)[1]*60))
print("------------------------\n")

print("HW 3.19")
print("------------------------")
print("S_W=%.12f ft^2\tV_MDV=%.12f ft/s\nP_in=%.12f hp\t(%.12f ft-lb/s)" % (HW3_19(0.0020482)[0],
                                                             HW3_19(0.0020482)[1],
                                                             HW3_19(0.0020482)[2],
                                                             HW3_19(0.0020482)[2]*550))
print("------------------------\n")

print("HW 3.20")
print("------------------------")
print("V_c=%.12f ft/s" % (HW3_20()))
print("------------------------\n")

print("HW 3.42")
print("------------------------")
print("50 mph tailwind")
print("V_BG=%.12f ft/s\t(%.12f mph)\nR_gmax=%.12f" % (HW3_42(-50*5280/3600)[0],
                                        HW3_42(-50*5280/3600)[0]*3600/5280,
                                        HW3_42(-50*5280/3600)[1]))
print("No wind")
print("V_BG=%.12f ft/s\t(%.12f mph)\nR_gmax=%.12f" % (HW3_42(0.)[0],
                                        HW3_42(0.)[0]*3600/5280,
                                        HW3_42(0.)[1]))
print("50 mph headwind")
print("V_BG=%.12f ft/s\t(%.12f mph)\nR_gmax=%.12f" % (HW3_42(50*5280/3600)[0],
                                        HW3_42(50*5280/3600)[0]*3600/5280,
                                        HW3_42(50*5280/3600)[1]))
print("------------------------\n")

print("HW 3.45")
print("------------------------")
HW3_45(7)
HW3_45(10)
HW3_45(17)
HW3_45(19)
print("------------------------\n")

print("HW 3.47")
print("------------------------")
print("Sea level")
HW3_47(0.0023769)
print("30,000 ft")
HW3_47(0.00089068)
print("------------------------\n")
