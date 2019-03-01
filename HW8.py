# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 13:06:23 2019

@author: Christian
"""

import numpy as np


def HW5_24(K_ell=0.992, K_gamma=0.797, eps_sbv=-0.035):
    # K_v = 1.0501, K_b = 0.7479, K_s = 1.0134, K_beta = 0.2834, CLw = 0.5230
    Sw = 950
    bw = 75
    CLwa = 4.3
    RTw = 0.35
    lamda_w = 15*np.pi/180.
    lw = 0.7
    hw = -3.8
    Sh = 207
    bh = 31
    CLha = 4.
    lh = 43.2
    hh = 13.7
    etah = 1.
    eps_dah = .3302
    Sv = 138
    cv = 6.9
    CLva = 3.4
    lv = 43.2
    hv = 7.
    etav = 1.
    RAw = bw*bw/Sw
    gamma_w = 4.2*np.pi/180.

    # Small angles in 5.6.13
    delCell_bgamw = -2.*K_ell*K_gamma*CLwa*gamma_w/(3.*np.pi)
    print("delCell_bgamw (0.2 pt):",delCell_bgamw)
    gamma_f = -3.5*np.pi/180.
    delCell_bf = -2.*K_ell*K_gamma*CLwa*gamma_f/(3.*np.pi)
    print("delCell_bf (0.2 pt):",delCell_bf)
    CLw = 0.5230
    delCell_blam = -2.*K_ell*CLw*np.tan(lamda_w)/(3.*np.pi)
    print("delCell_blam (0.2 pt):",delCell_blam)
    delCell_bv = -etav*Sv*hv*(1. - eps_sbv)*CLva/(Sw*bw)
    print("delCell_bv (0.2 pt):",delCell_bv)
    delCell_bh = -0.08*etav*Sv*bh*(1. - eps_sbv)*CLva/(Sw*bw)
    print("delCell_bh (0.2 pt):",delCell_bh)
    Cell_b = delCell_bgamw + delCell_bf + delCell_blam + delCell_bv + delCell_bh
    print("Cell_b (1 pt):", Cell_b)



def HW5_26():
    Sw = 327.5
    bw = 52
    CLa = 5.01
    V_inf = 352
    W = 15480
    dp = 138./12.
    hp = 0.
    omega_2pi = 1350/60
    eps_dap = -0.175
    Cnpa = 0.046
    CD = 0.0466
    Clda = -0.131
    Cnda = 0.013
    Cldr = 0.006
    Cndr = -0.052
    V_max = 400.*5280/3600
    rho_sl = 0.0023769
    rho_25k = 0.0010663
    etap = 0.8
    CL_400 = W/(0.5*rho_sl*V_max*V_max*Sw)
    CL_240 = W/(0.5*rho_25k*V_inf*V_inf*Sw)
    print("CL_400:", CL_400)
    print("CL_240:", CL_240)
    a0p = (CL_240 - CL_400)/CLa
    print("a0p:", a0p)
    eps_d0p = CL_240*eps_dap/CLa
    print("eps_d0p:", eps_d0p)
    T = 0.5*rho_25k*V_inf*V_inf*Sw*CD
    print("Thrust:", T)
    y_bp = 8.
    deln_p = -y_bp*T + rho_25k*omega_2pi*omega_2pi*dp*dp*dp*dp*dp*Cnpa*(a0p - eps_d0p)
    delCn_p = deln_p/(0.5*rho_25k*V_inf*V_inf*Sw*bw)
    print("delCn_p (0.5 pts):", delCn_p)
    dell_p = -T*V_inf/(etap*2.*np.pi*-omega_2pi)
    delCell_p = dell_p/(0.5*rho_25k*V_inf*V_inf*Sw*bw)
    print("delCell_p (0.5 pts):", delCell_p)
    A = [[Clda, Cldr],
         [Cnda, Cndr]]
    b = [-delCell_p, -delCn_p]
    x = np.linalg.solve(A, b)
    da = x[0]
    dr = x[1]
    print("da (0.5 pts):", da*180./np.pi, "deg\n\t     ", da, "rad")
    print("dr (0.5 pts):", dr*180./np.pi, "deg\n\t     ", dr, "rad")


def HW5_27():
    Sw = 327.5
    bw = 52.
    CLa = 5.01
    V_inf = 200
    W = 15480
    Pb = 1425*550
    etap = 0.7
    dp = 138/12
    hp = 0.
    omega_2pi = -1500/60
    Cnpa = 0.0089
    eps_dap = -0.175
    CYb = -0.512
    CYda = 0.
    CYdr = 0.115
    Clb = -0.025
    Clda = -0.131
    Cldr = 0.006
    Cnb = 0.083
    Cnda = 0.013
    Cndr = -0.052
    rho_sl = 0.0023769
    V_max = 400*5280/3600
    CL = W/(0.5*rho_sl*V_inf*V_inf*Sw)
    CL_400 = W/(0.5*rho_sl*V_max*V_max*Sw)
    a0p = (CL - CL_400)/CLa
    print("(0.5 pts. for these three components)")
    print("a0p:", a0p)
    eps_d0p = CL*eps_dap/CLa
    print("eps_d0p:", eps_d0p)
    T = etap*Pb/V_inf
    print("T:", T)
    ybp = 8.
    deln_p = ybp*T + rho_sl*omega_2pi*omega_2pi*dp*dp*dp*dp*dp*Cnpa*(a0p - eps_d0p)
    delCn_p = deln_p/(0.5*rho_sl*V_inf*V_inf*Sw*bw)
    print("delCn_p (0.5 pts):", delCn_p)
    dell_p = -Pb/(2.*np.pi*omega_2pi)
    delCell_p = dell_p/(0.5*rho_sl*V_inf*V_inf*Sw*bw)
    print("delCell_p (0.5 pts):", delCell_p)
    A = [[CYb, CYda, CYdr],
         [Clb, Clda, Cldr],
         [Cnb, Cnda, Cndr]]
    delCYp = 0.
    bank = 0.
    CW = CL
    b = [-delCYp - CW*bank, -delCell_p, -delCn_p]
    x = np.linalg.solve(A, b)
    beta = x[0]
    da = x[1]
    dr = x[2]
    print("beta (0.5 pts):", beta*180./np.pi, "deg\n\t       ", beta, "rad")
    print("da (0.5 pts):", da*180./np.pi, "deg\n\t     ", da, "rad")
    print("dr (0.5 pts):", dr*180./np.pi, "deg\n\t     ", dr, "rad")


def V_bodyfixed(V, a, b, experiment=False):
    SA = np.sin(a)
    CA = np.cos(a)
    SB = np.sin(b)
    CB = np.cos(b)
    if not experiment:
        denom = np.sqrt(1. - SA*SA*SB*SB)
        V_xb = V*CA*CB/denom
        V_yb = V*CA*SB/denom
        V_zb = V*SA*CB/denom
        return V_xb, V_yb, V_zb
    else:
        V_xb = V*CA*CB
        V_yb = V*SB
        V_zb = V*SA*CB
        return V_xb, V_yb, V_zb


def HW7_3():
    V = 400
    alpha = 20*np.pi/180.
    b_e = 10*np.pi/180.
    V_xb, V_yb, V_zb = V_bodyfixed(V, alpha, b_e, experiment=True)
    print("(0.5 pts. for these three components)")
    print("V_xb:", V_xb)
    print("V_yb:", V_yb)
    print("V_zb:", V_zb)



def HW7_4():
    V = 800
    alpha = 75*np.pi/180.
    b_e = 60*np.pi/180.
    V_xb, V_yb, V_zb = V_bodyfixed(V, alpha, b_e)
    print("(0.5 pts. for these three components)")
    print("V_xb:", V_xb)
    print("V_yb:", V_yb)
    print("V_zb:", V_zb)

def HW7_5():
    V = 800
    alpha = 75*np.pi/180.
    b_e = 60*np.pi/180.
    V_xb, V_yb, V_zb = V_bodyfixed(V, alpha, b_e, experiment=True)
    print("(0.5 pts. for these three components)")
    print("V_xb:", V_xb)
    print("V_yb:", V_yb)
    print("V_zb:", V_zb)


print("5.24 (2.0 pts)")
print("-------------------------")
HW5_24()
print("\n5.26 (2.0 pts)")
print("-------------------------")
HW5_26()
print("\n5.27 (3.0 pts)")
print("-------------------------")
HW5_27()
print("\n7.3 (1.0 pt)")
print("-------------------------")
HW7_3()
print("\n7.4 (1.0 pt)")
print("-------------------------")
HW7_4()
print("\n7.5 (1.0 pt)")
print("-------------------------")
HW7_5()

