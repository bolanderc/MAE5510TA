#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 10:30:42 2019

@author: christian
"""

import json
import numpy as np


def HW1_37(V):
    c = 5.5
    RA = 6
    b = RA*c
    print("Semispan:", b/2, "ft")
    file = json.load(open("HW1_37Wing.json"))
    CL = file['total']['MyAirplane']['CL']
    rho = 0.0023769
    L = 0.5*rho*V*V*b*c*CL
    print("Lift:", L, "lbf")


def HW1_38(V):
    RA = 4.0
    S = 36
    b = np.sqrt(RA*S)
    print("Semispan:", b/2., "ft")
    RT = 0.5
    # Root chord found by 0.5(cr + ct)b = S
    c_r = (2.*S)/(b*(1. + RT))
    c_t = RT*c_r
    print("Root chord:", c_r, "ft")
    print("Tip chord:", c_t, "ft")
    c = (c_r + c_t)/2.
    # Sweep at the quarter chord
    sweep = np.arctan2(0.75*(c_r - c_t), b/2)*180/np.pi
    print("Quarter Chord Sweep:", sweep, "deg", "(", sweep*np.pi/180., " rad)")
    file = json.load(open("HW1_38Wing.json"))
    CL = file['total']['MyAirplane']['CL']
    rho = 0.0023769
    L = 0.5*rho*V*V*b*c*CL
    print("Lift:", L, "lbf")


def HW1_39(V):
    file39 = json.load(open("HW1_39WingTail.json"))
    CL_wing = file39['total']['Wing_3_right']['CL']*2.
    rho = 0.0023769
    S = 181.5
    L_wing = 0.5*rho*V*V*S*CL_wing
    print("Wing Lift:", L_wing, "lbf")
    CL_tail = file39['total']['Wing_2_right']['CL']*2.
    L_tail = 0.5*rho*V*V*S*CL_tail
    print("Tail Lift:", L_tail, "lbf")


def HW4_1():
    W = 2000
    rho = 0.0014962
    V_inf = 220
    """
    Unknowns are Sw, lw, a_o_w, St, lt, a_o_t
    (1) The lift must support the weight, L_w + L_t = W
    (2) The pitching moment about the CG must be zero, Lwlw + Ltlt = 0
    (3) The airplane’s neutral point is to be 0.5 ft aft of the CG.
        Thus from Eq. (4.4.3), Swlw + Stlt = 0.5(Sw + St)
    3 eqs. 6 unknowns
    Constraint on fuselage length of 15 ft.
    - Let's design for minimum drag with no lift on the tail.
    - Mount the tail 7.5 ft aft of the CG
    - Assume value for e, CDo, and RAw
    Use CL definition to find Sw
    No downwash means that mounting angle for tail must be zero
    Use equations 1-3 for mounting of main wing, location of main wing, and St

    """
    l_t = 7.5
    print("l_t:", l_t, "ft")
    e = 0.8
    CDo = 0.02
    RA_w = 6.
    Sw = W/(np.sqrt(np.pi*e*CDo*RA_w)*rho*V_inf*V_inf*0.5)
    print("S_w:", Sw, "ft^2")
    a_o_t = 0.
    print("a_o_t:", a_o_t, "deg")
    # From (1)
    a_o_w = (W/(rho*V_inf*V_inf*np.pi))/Sw
    print("a_o_w:", a_o_w*180/np.pi, "deg", "(", a_o_w, " rad)")
    # From (2)
    l_w = 0.
    print("l_w:", l_w, "ft")
    # From (3)
    St = 0.5*Sw/(l_t - 0.5)
    print("S_t:", St, "ft^2")


def HW4_3(V):
    V = V*5280/3600  # mph to fps
    rho = 0.0023769
    Sw = 244
    bw = 36.83
    CLw_a = 4.62
    aLO_w = -2.2*np.pi/180.
    Cmw = -0.053
    Sh = 31
    CLh_a = 4.06
    lh_lw = 18.16
    W = 8375
    # With no elevator deflection, the pitching moment for the symmetric
    # horizontal tail is zero and with the CG at
    # the wing quarter-chord, l_w is also zero.
    lw = 0
    lh = lh_lw
    # Assume tail efficiency of 1.0
    eta_h = 1.0
    cw = Sw/bw
    a0_h = Sw*cw*Cmw/(Sh*lh*eta_h*CLh_a)
    print("a0_h:", a0_h, "(", a0_h*180./np.pi, " deg)")
    # Lift on wing and horizontal tail must support weight
    a0_w = ((W/(0.5*rho*V*V*Sw*CLw_a)) - (Sh*eta_h*CLh_a)*a0_h/(Sw*CLw_a) +
            aLO_w)
    print("a0_w:", a0_w, "(", a0_w*180./np.pi, " deg)")
    # Eq. (4.3.30)
    eps_da = 0.
    dCm_da = (-lw/cw)*CLw_a - (Sh*lh)/(Sw*cw)*eta_h*CLh_a*(1 - eps_da)
    print("dCm/da:", dCm_da)


def HW4_4():
    Sw = 244
    CLw_a = 4.62
    Sh = 31
    CLh_a = 4.06
    eta_t = 1.0
    eps_da = 0.0
    dCL_da = CLw_a + (Sh/Sw)*eta_t*CLh_a*(1 - eps_da)
    dCm_da = -1.414
    # Eq. (4.4.13)
    stat_marg = -(dCm_da/dCL_da)
    print("l_np/cw:", stat_marg)


print("1.37")
print("-------------------------")
HW1_37(176)
print("\n1.38")
print("-------------------------")
HW1_38(176)
print("\n1.39")
print("-------------------------")
HW1_39(176)
print("\n4.1")
print("-------------------------")
HW4_1()
print("\n4.3")
print("-------------------------")
HW4_3(200)
print("\n4.4")
print("-------------------------")
HW4_4()
