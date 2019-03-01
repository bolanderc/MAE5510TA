#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 15:23:19 2019

@author: christian
"""

import numpy as np
import scipy.optimize as optimize


def K_s(xbar, ybar, sweep, Kb):
    rbar = np.sqrt(xbar*xbar + ybar*ybar)
    sbar = Kb*np.tan(sweep)
    tbar = np.sqrt((xbar - sbar)*(xbar - sbar) + ybar*ybar + Kb*Kb)
    t0bar = np.sqrt(xbar*xbar + ybar*ybar + Kb*Kb)
    C1 = (xbar - sbar)/(tbar)
    C2 = (xbar*(rbar + tbar)*(t0bar*t0bar - xbar*xbar))/(rbar*tbar*(rbar*tbar + rbar*rbar - xbar*sbar))
    num = 1. + C1 + C2
    denom = 1. + (xbar*(rbar*rbar + t0bar*t0bar - xbar*xbar))/(rbar*rbar*t0bar)
    return num/denom


def K_p(xbar, ybar, Kb):
    S1 = (2.*Kb*Kb)/(np.pi*np.pi*(ybar*ybar + Kb*Kb))
    C1 = xbar*(xbar*xbar + 2.*ybar*ybar + Kb*Kb)
    C2 = (xbar*xbar + ybar*ybar)*np.sqrt(xbar*xbar + ybar*ybar + Kb*Kb)
    return S1*(1. + C1/C2)


def HW4_7():
    Sw = 244  # ft^2
    bw = 36.83  # ft
    CLwa = 4.62
    lw = 0.  # ft
    Sh = 31.  # ft^2
    bh = 10.64  # ft
    CLha = 4.06
    lh = 18.16  # ft
    etah = 1.
    cw = Sw/bw
    # In problem 4.5
    RAw = bw*bw/Sw
    xbar = lh/(bw/2.)
    ybar = 0.
    # Fig. 4.5.2 and 4.5.3 yeild
    Kv = 1.0
    Kb = np.pi/4.
    # No sweep
    Ks = 1.
    # Eq. 4.5.6
    Kp = K_p(xbar, ybar, Kb)
    epsda = (Kv*Kp*Ks*CLwa)/(Kb*RAw)
    print(u"\u03B5_d,a :", epsda)
    CLa = CLwa + (Sh/Sw)*etah*CLha*(1 - epsda)
    print('CLa : ', CLa)
    Cma = -(lw/cw)*CLwa - (Sh*lh)/(Sw*cw)*etah*CLha*(1 - epsda)
    print('Cma : ', Cma)
    stat_marg = -Cma/CLa
    print("Static Margin:", stat_marg)


def HW4_12(V):
    V *= (5280/3600)
    Sw = 244  # ft^2
    bw = 36.83  # ft
    CLwa = 4.62
    Cmw = -0.053
    rho = 0.0023769
    aL0w = -2.2  # deg
    a0w = 2.0  # deg
    aL0w *= (np.pi/180)  # rad
    a0w *= (np.pi/180)  # rad
    lw = 0.  # ft
    hw = -2.6  # ft
    CD0w = 0.01
    CD0Lw = 0.
    ew = 1.
    cw = Sw/bw
    RAw = bw*bw/Sw
    Sh = 31.  # ft^2
    bh = 10.64  # ft
    CLha = 4.06
    Cmh0 = 0.
    aL0h = 0.
    a0h = 0.
    lh = 18.16  # ft
    hh = 0.
    CD0h = 0.01
    CD0Lh = 0.
    eh = 1.
    W = 8375
    epse = 0.6
    Cmhde = -0.55
    lp = -9.33  # ft
    hp = 0.
    RAh = bh*bh/Sh
    etah = 1.0
    ch = Sh/bh
    # Using method in 4.5 downwash is estimated
    xbar = (lh)/(bw/2)
    ybar = (hh - hw)/(bw/2)
    Kv = 1.0
    Kb = 0.7864
    Kp = K_p(xbar, ybar, Kb)
    Ks = 1.0

    def force_moment(c):
        eps_d = (Kv*Kp*Ks*CLwa)*(c[0] + a0w - aL0w)/(Kb*RAw)
        CLw = CLwa*(c[0] + a0w - aL0w)
        CLh = CLha*(c[0] + a0h - aL0h - eps_d + epse*c[1])
        CDw = CD0w + CD0Lw*CLw + (CLw*CLw)/(np.pi*ew*RAw)
        CDh = CD0h + CD0Lh*CLh + (CLh*CLh)/(np.pi*eh*RAh)
        Cmh = Cmh0 + Cmhde*c[1]
        err1 = (CLw*np.cos(c[0]) + CDw*np.sin(c[0]) +
                (Sh/Sw)*(CLh*np.cos(c[0] - eps_d) + CDh*np.sin(c[0] - eps_d)) -
                (W*np.cos(c[0]))/(.5*rho*V*V*Sw))
        err2 = (Cmw - (lw/cw)*(CLw*np.cos(c[0]) + CDw*np.sin(c[0])) +
                (hw/cw)*(CDw*np.cos(c[0]) - CLw*np.sin(c[0])) +
                (Sh*ch)/(Sw*cw)*Cmh -
                (Sh*lh)/(Sw*cw)*(CLh*np.cos(c[0] - eps_d) +
                                 CDh*np.sin(c[0] - eps_d)))
        error = (abs(err1) + abs(err2))/2.
        return error
    res = optimize.minimize(force_moment, [0.5, -0.5]).x
    alpha = res[0]
    d_e = res[1]
    print("Velocity : ", round(V*(3600/5280), 2), "mph")
    print("alpha :", np.rad2deg(alpha))
    print("d_e :", np.rad2deg(d_e))
    eps_d = (Kv*Kp*Ks*CLwa)*(alpha + a0w - aL0w)/(Kb*RAw)
    Lh = 0.5*rho*V*V*Sh*etah*CLha*(alpha + a0h - eps_d + epse*d_e)
    print("Lift on tail :", Lh)


def HW4_23(V):
    V *= (5280/3600)
    rho = 8.9068e-4
    Sw = 950
    bw = 75
    RAw = bw*bw/Sw
    cw = Sw/bw
    CLwa = 4.3
    aL0w = -2.0*np.pi/180.
    Cmw = -0.05
    Sh = 219
    bh = 36
    RAh = bh*bh/Sh
    ch = Sh/bh
    CLha = 4.3
    aL0h = -2.0*np.pi/180.
    Cmh0 = -0.05
    lw_lh = 42.5
    eps_e = 0.51
    Cmhde = -0.49
    sweep = 15*np.pi/180.
    stat_mar = 0.05
    W = 73000
    xbar = -(lw_lh)/(bw/2.)
    ybar = 0.
    Kv = 1.0501
    Kb = 0.7479
    Kp = K_p(xbar, ybar, Kb)
    Ks = K_s(xbar, ybar, sweep, Kb)
    epsda = (Kv*Kp*Ks*CLwa)/(Kb*RAw)
    print("eps_da :", epsda)
    CLa = CLwa + (Sh/Sw)*CLha*(1. - epsda)
    Cma = -stat_mar*CLa
    lw = (-Cma*cw + (Sh/Sw)*CLha*(1. - epsda)*lw_lh)/(CLwa + (Sh/Sw)*CLha*(1. - epsda))
    lh = lw - lw_lh
    A = np.zeros((2, 2))
    A[0, 0] = 1.
#    A[0, 0] = 1.
    A[1, 0] = 1.
#    A[1, 0] = -(lw/cw)
    A[0, 1] = (Sh*CLha)/(Sw*CLwa - epsda*Sh*CLha)
#    A[0, 1] = (Sh/Sw)
    A[1, 1] = (Sh*lh*CLha)/(Sw*lw*CLwa - epsda*Sh*lh*CLha)
#    A[1, 1] = -Sh*lh/(Sw*cw)
    B = np.zeros(2)
    B[0] = ((W/(0.5*rho*V*V)) + Sh*CLha*aL0h)/(Sw*CLwa - epsda*Sh*CLha) + aL0w
#    B[0] = W/(0.5*rho*V*V*Sw)
    B[1] = (Sw*cw*Cmw + Sh*ch*Cmh0 + Sh*lh*CLha*aL0h)/(Sw*lw*CLwa - epsda*Sh*lh*CLha) + aL0w
#    B[1] = -Cmw - (Sh*ch*Cmh0)/(Sw*cw)
    mounts = np.linalg.solve(A, B)
#    mounts[0] = mounts[0]/CLwa + aL0w
#    mounts[1] = (mounts[1] + CLha*(aL0h + epsda*aL0h))/(CLha*(1 + epsda))
    print("CL_a :", CLa)
    print("Cm_a :", Cma)
    print("l_w :", lw)
    print("l_h :", lh)
    print("alpha_0w :", np.rad2deg(mounts[0]), "deg (", mounts[0], "rad)")
    print("alpha_0h :", np.rad2deg(mounts[1]), "deg (", mounts[1], "rad)")


print("4.7")
print("-------------------------")
HW4_7()
print("\n4.12")
print("-------------------------")
HW4_12(100)
print("-------------------------")
HW4_12(200)
print("-------------------------")
HW4_12(300)
print("\n4.23")
print("-------------------------")
HW4_23(400)

print(K_p(-0.5067, 0., np.pi/4.))

#4.7 (/3) :
#4.12 (/3.5) :
#4.23 (/3.5) :
#Total (/10) :
