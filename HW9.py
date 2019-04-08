#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 14:57:10 2019

@author: christian
"""

import numpy as np
import scipy.optimize as optimize


def matprint(mat, fmt="g"):
    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
        print("")


def HW7_32():
    ct_crw = 0.286
    Sv = 875
    xbv = -96
    zbv = -23
    CLvB = 3.3
    CYB = -0.890
    ClB = -0.144
    CnB = 0.182
    Sw = 5500
    bw = 196
    cw = 28
    W = 636600
    CLa = 5.5
    CLwa = 4.67
    Sh = 1300
    xbh = -100
    lwt = 55
    CLha = 3.5
    Cma = -1.26
    rho = 0.00089068
    Vo = 660
    M = 0.663
    Fybv = rho*Vo*Sw*CYB/2
    print("Fybv (0.33 pts):", Fybv)
    Mxbv = rho*Vo*Sw*bw*ClB/2.
    print("Mxbv (0.33 pts):", Mxbv)
    Mzbv = rho*Vo*Sw*bw*CnB/2.
    print("Mzbv (0.33 pts):", Mzbv)
    # Assuming all lift on the wing
    CLw = W/(0.5*rho*Vo*Vo*Sw)
    RAw = bw*bw/Sw
    K_lp = 0.79
    Fybp = 0.
    print("Fybp (0.33 pts):", Fybp)
    Mxbp = -rho*Vo*Sw*bw*bw*K_lp*CLwa/32.
    print("Mxbp (0.33 pts):", Mxbp)
    Mzbp = -rho*Vo*Sw*bw*bw*(1. - 3.*K_lp*CLwa/(np.pi*RAw))*CLw/32.
    print("Mzbp (0.33 pts):", Mzbp)
    etav = 1.
    Fybr = -etav*(0.5*rho*Vo*Sv)*xbv*CLvB
    print("Fybr (0.33 pts):", Fybr)
    Mxbr = (0.5*rho*Vo)*((1. - 3.*K_lp*CLwa/(np.pi*RAw))*Sw*bw*bw*CLw/8. +
                         etav*Sv*xbv*zbv*CLvB)
    print("Mxbr (0.33 pts):", Mxbr)
    Mzbr = -etav*(0.5*rho*Vo*Sv)*xbv*xbv*CLvB
    print("Mzbr (0.33 pts):", Mzbr)


def HW7_34():
    Sw = 5500
    bw = 196
    cw = 28
    xbw = 0.
    W = 636600
    CLa = 5.5
    CLwa = 4.67
    Sh = 1300
    xbh = -100
    lwt = 55
    CLha = 3.5
    Cma = -1.26
    rho = 0.00089068
    Vo = 660
    M = 0.663
    CLw = 0.597
    RAw = 6.98
    Klp = 0.79
    RTw = 0.286
    Sv = 875
    xbv = -96
    zbv = -23
    CLvB = 3.3
    CYB = -0.89
    ClB = -0.144
    CnB = 0.182
    Ixxb = 1.836e7
    Izzb = 4.954e7
    Ixzb = 7.33e5
    etav = 1.
    g = 32.2
    theta_o = 0.
    CYp = 0
    Clp = -Klp*CLwa/8
    Cnp = -(1. - 3.*Klp*CLwa/(np.pi*RAw))*CLw/8.
    CYr = -etav*2.*Sv*xbv*CLvB/(Sw*bw)
    Clr = (1. - 3.*Klp*CLwa/(np.pi*RAw))*CLw/4. + etav*2.*Sv*xbv*zbv*CLvB/(Sw*bw*bw)
    Cnr = -etav*2.*Sv*xbv*xbv*CLvB/(Sw*bw*bw)
    ixz = Ixzb/Ixxb
    izx = Ixzb/Izzb
    Rgy = g*bw/(2.*Vo*Vo)
    Cm = rho*Sw*bw/(4.*W/g)
    Cx = rho*Sw*bw*bw*bw/(8.*Ixxb)
    Cz = rho*Sw*bw*bw*bw/(8.*Izzb)
    RyB = Cm*CYB
    RlB = Cx*ClB
    RnB = Cz*CnB
    Ryp = Cm*CYp
    Rlp = Cx*Clp
    Rnp = Cz*Cnp
    Ryr = Cm*CYr
    Rlr = Cx*Clr
    Rnr = Cz*Cnr
    LHS = np.zeros((6, 6))
    np.fill_diagonal(LHS, 1.)
    LHS[2, 1] = -izx
    LHS[1, 2] = -ixz
    RHS = np.zeros((6, 6))
    RHS[0, 0] = RyB
    RHS[0, 1] = Ryp
    RHS[0, 2] = Ryr - 1.
    RHS[0, 4] = Rgy*np.cos(theta_o)
    RHS[1, 0] = RlB
    RHS[1, 1] = Rlp
    RHS[1, 2] = Rlr
    RHS[2, 0] = RnB
    RHS[2, 1] = Rnp
    RHS[2, 2] = Rnr
    RHS[3, 0] = 1.
    RHS[3, 5] = np.cos(theta_o)
    RHS[4, 1] = 1.
    RHS[4, 2] = np.tan(theta_o)
    RHS[5, 2] = 1./np.cos(theta_o)
    print("LHS (1 point):")
    print("-----------------")
    matprint(LHS)
    print("ixz:", ixz)
    print("izx:", izx)
    print("\nRHS (3 points):")
    print("-----------------")
    matprint(RHS)
    print("RyB:", RyB)
    print("Ryp:", Ryp)
    print("Ryr:", Ryr)
    print("Rgy:", Rgy)
    print("RlB:", RlB)
    print("Rlp:", Rlp)
    print("Rlr:", Rlr)
    print("RnB:", RnB)
    print("Rnp:", Rnp)
    print("Rnr:", Rnr)


def HW7_35():
    Ixxb = 1.888e7
    Iyyb = 3.31e7
    Izzb = 4.902e7
    Ixzb = 4.075e6
    def Eq7_8_8(phi):
        err = Ixzb*(np.cos(phi)**2 - np.sin(phi)**2) + (Izzb - Ixxb)*np.cos(phi)*np.sin(phi)
        return np.abs(err)
    x = optimize.minimize(Eq7_8_8, 0.).x
    print("phi_p:", x*180/np.pi, "deg\n      ", x, "rad")


print("\n7.32 (3.0 pts)")
print("-------------------------")
HW7_32()
print("\n7.34 (4.0 pts)")
print("-------------------------")
HW7_34()
print("\n7.35 (3.0 pts)")
print("-------------------------")
HW7_35()

#7.32 (/3) :
#7.34 (/4) :
#7.35 (/3) :
#Total (/10) :
