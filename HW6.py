#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 09:52:53 2019

@author: christian
"""

import numpy as np


def K_p(xbar, ybar, Kb):
    S1 = (2.*Kb*Kb)/(np.pi*np.pi*(ybar*ybar + Kb*Kb))
    C1 = xbar*(xbar*xbar + 2.*ybar*ybar + Kb*Kb)
    C2 = (xbar*xbar + ybar*ybar)*np.sqrt(xbar*xbar + ybar*ybar + Kb*Kb)
    return S1*(1. + C1/C2)


def HW4_10():
    Sw = 244
    bw = 36.83
    CLwa = 4.62
    Cmw = -0.053
    aL0w = -2.20*np.pi/180.
    a0w = 2.0*np.pi/180.
    lw = 0.
    Sh = 31
    bh = 10.64
    CLha = 4.06
    Cmh0 = 0.0
    aL0h = 0.0
    a0h = 0.0
    lh = 18.16
    W = 8375
    eps_e = 0.6
    Cmhde = -0.55
    etah = 1.
    eps_d0 = 0.0358
    eps_da = 0.4883
    Sf = 15
    cf = 32.67
    l_CG = 11
    l_maxCS = 16

    cw = Sw/bw
    df = 2.0*np.sqrt(Sf/np.pi)
    lf = (l_maxCS/2.) - l_CG
    print("eps_d0:", eps_d0)
    print("eps_da:", eps_da)
    print("df:", df)
    print("lf:", lf)
    CLa_wh = CLwa + (Sh/Sw)*etah*CLha*(1. - eps_da)
    print("CLa_wh:", CLa_wh)
    Cma_wh = (-lw/cw)*CLwa - (Sh*lh/(Sw*cw))*etah*CLha*(1.0 - eps_da)
    print("Cma_wh:", Cma_wh)
    delCma_f = -2.*(Sf*lf/(Sw*cw))*(1.0 - 1.76*(df/cf)**1.5)
    print("delCma_f:", delCma_f)
    stat_marg_approx = -(Cma_wh + delCma_f)/(CLa_wh)
    print("static margin close approx using Eq. (4.4.13):", stat_marg_approx)
    stat_marg = -(Cma_wh + delCma_f)/(CLa_wh - delCma_f*(cw/lf))
    print("static margin following Ex. 4.9.1:", stat_marg)


def HW4_11():
    Sw = 244
    bw = 36.83
    CLwa = 4.62
    Cmw = -0.053
    aL0w = -2.20*np.pi/180.
    a0w = 2.0*np.pi/180.
    lw = 0.
    Sh = 31
    bh = 10.64
    CLha = 4.06
    Cmh0 = 0.0
    aL0h = 0.0
    a0h = 0.0
    lh = 18.16
    W = 8375
    eps_e = 0.6
    Cmhde = -0.55
    etah = 1.
    eps_d0 = 0.0358
    eps_da = 0.4883
    Sf = 15
    cf = 32.67
    l_CG = 11
    l_maxCS = 16


    cw = Sw/bw
    df = 2.0*np.sqrt(Sf/np.pi)
    lf = (l_maxCS/2.) - l_CG
    dp = 132./12.  # ft
    omega = 1300./60.  # rps
    lp = -9.33
    V_inf = 200*5280/3600  # ft/s
    CNp_a = 0.27

    CLa_wh = CLwa + (Sh/Sw)*etah*CLha*(1. - eps_da)
    Cma_wh = (-lw/cw)*CLwa - (Sh*lh/(Sw*cw))*etah*CLha*(1.0 - eps_da)
    delCma_f = -2.*(Sf*lf/(Sw*cw))*(1.0 - 1.76*(df/cf)**1.5)
    J = V_inf/(omega*dp)
    print("Advance Ratio:", J)

    RAw = (bw*bw/Sw)
    xbar = (lp - lw)/(bw/2.)
    ybar = 0.0
    Kv = 1.0
    Ks = 1.0
    Kb = np.pi/4.
    Kp = K_p(xbar, ybar, Kb)
    eps_da_p = (Kv*Kp*Ks*CLwa/(Kb*RAw))
    print("Kp:", Kp)
    print("eps_da_p:", eps_da_p)
    delCma_p = -(2.*dp*dp*lp*(1. - eps_da_p)*CNp_a)/(Sw*cw*J*J)
    print("delCma_p:", delCma_p)
    stat_marg_approx = -(Cma_wh + delCma_f + delCma_p)/(CLa_wh)
    print("static margin close approx using Eq. (4.4.13):", stat_marg_approx)
    stat_marg = -(Cma_wh + delCma_f + delCma_p)/(CLa_wh - delCma_f*(cw/lf) -
                                                 delCma_p*(cw/lp))
    print("static margin following Ex. 4.9.1:", stat_marg)


print("4.10")
print("-------------------------")
HW4_10()
print("\n4.11")
print("-------------------------")
HW4_11()


#Ex. 4.8.2 (/3.5):
#4.10 (/3):
#4.11 (/3.5):
#Total (/10):
