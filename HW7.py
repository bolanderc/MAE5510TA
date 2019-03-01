#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 12:05:36 2019

@author: christian
"""

import numpy as np


def HW5_3():
    Sw = 244
    bw = 36.83
    cw = Sw/bw
    CLwa = 4.62
    Cmw = -0.053
    aL0w = -2.2*np.pi/180.
    a0w = 2.0*np.pi/180.
    lw = 0.
    Sh = 31.
    bh = 10.64
    ch = Sh/bh
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
    df = 2.*np.sqrt(Sf/np.pi)
    cf = 32.67
    lcp = 16.
    lcg = 11.
    lf = lcp/2. - lcg
    dp = 132./12.
    omega = 1300./60.
    lp = -9.33
    V_inf = 200*5280/3600
    CNpa = 0.27
    Cn_b = 0.1
    CLva = 3.5
    """
    Assumptions
    """
    eps_da_p = 0.
    eps_sb_v = 0.
    etav = 1.
    """
    Solution
    """
    delCma_f = -2.*(Sf*lf/(Sw*cw))*(1.0 - 1.76*(df/cf)**1.5)
    delCn_B_f = -delCma_f*cw/bw
    print("delCn_B_f (2pts):", delCn_B_f)
    J = V_inf/(omega*dp)
    print("J (2pts):", J)
    delCn_B_p = 2.*dp*dp*lp*(1. - eps_da_p)*CNpa/(Sw*bw*J*J)
    print("delCn_B_p (2pts):", delCn_B_p)
    delCn_B_v = Cn_b - delCn_B_f - delCn_B_p
    print("delCn_B_v (2pts):", delCn_B_v)
    VTVR = delCn_B_v/(etav*CLva*(1. - eps_sb_v))
    print("VTVR (2pts):", VTVR)


def HW5_4():
    Sv = 21.
    CLva = 3.5
    lv = 18.16
    etav = 1.
    eps_sBv = 0.
    cv = 3.8
    eps_r = 0.6
    Cmvdr = -0.58
    Sw = 244.
    bw = 36.83
    CLwa = 4.62
    Cmw = -0.053
    W = 8375
    aL0w = -2.2*np.pi/180.
    a0w = 2.*np.pi/180.
    lw = 0.
    Sh = 31.
    bh = 10.64
    CLha = 4.06
    Cmh0 = 0.
    aL0h = 0.
    a0h = 0.
    lh = 18.16
    dp = 132./12.
    omega = 1400./60.
    lp = -9.33
    V_inf = 55.*5280./3600.
    eps_dap = -0.1809
    Cnpa_mag = 0.0785
    a_prop = 11
    Cndr = -etav*Sv*lv*(eps_r*CLva - (cv/lv)*Cmvdr)/(Sw*bw)
    print("Cndr (1 pt):", Cndr)
    J = V_inf/(omega*dp)
    print("J (1 pt):", J)
    Cna_mag = 2.*dp*dp*dp*(1. - eps_dap)*Cnpa_mag/(Sw*bw*J*J)
    print("Cna_mag (1 pt):", Cna_mag)
    dr_mag = np.abs(Cna_mag*a_prop/Cndr)
    print("dr_mag (2 pts):", dr_mag)



print("5.3 (10 pts)")
print("[One point for equation, one for answer]")
print("-------------------------")
HW5_3()
print("\n5.4 (Extra Credit: 5 pts)")
print("[Half point for equation, remainder for answer]")
print("-------------------------")
HW5_4()
