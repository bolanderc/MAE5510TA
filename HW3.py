#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 10:21:15 2019

@author: christian

Christian Bolander
MAE 5510
HW 3 Problem 3.55
"""
import numpy as np
import matplotlib.pyplot as plt


def CD(CDo, CDoL, Cl, hw, bw, e, RA):
    C1 = CDo
    C2 = CDoL*Cl
    C3n = (Cl**2)*((16*hw/bw)**2)
    C3d = (1+((16*hw/bw)**2))*(np.pi*e*RA)
    return C1 + C2 + (C3n/C3d)


def T(V):
    return 199000 - (177.6*V) + (0.1193*(V**2))


def VLO(Clmax, W, Sw, p):
    C1 = np.sqrt(2/Clmax)
    C2 = np.sqrt((W/Sw)/p)
    return 1.1*C1*C2


def VOC(Clmax, W, Sw, p):
    C1 = np.sqrt(2/Clmax)
    C2 = np.sqrt((W/Sw)/p)
    return 1.2*C1*C2


def K0(Ts, W, u):
    C1 = Ts/W
    return C1-u


def K1(Tp, W):
    C1 = Tp/W
    return C1


def K2(Tpp, W, p, Sw, Cl, u, Cd):
    C1 = Tpp/W
    C2n = p*(Cl*u - Cd)
    C2d = (2*W)/Sw
    return C1 + (C2n/C2d)


def fip1(K0, K1, VLO, K2):
    C1 = K0
    C2 = K1*VLO
    C3 = K2*(VLO**2)
    return C1+C2+C3


def fi(K0, K1, Vi, K2):
    C1 = K0
    C2 = K1*Vi
    C3 = K2*(Vi**2)
    return C1 + C2 + C3


def fpi(K1, K2, Vi):
    C1 = K1
    C2 = 2*K2*Vi
    return C1 + C2


def KT(Vi, K0, K1, K2, fip1, fi, Vip1, Kw):
    if K2 == 0.0:
        if K1 == 0.0:
            Kt = ((Vip1**2) - (Vi**2))/(2*K0)
        else:
            Kt = (K0/(K1**2))*np.log(fi/fip1) + ((Vip1-Vi)/K1)
    else:
        Kt = (1/(2*K2))*np.log(fip1/fi) - ((K1*Kw)/(2*K2))
    return Kt


def KW(Kr, fpip1, fpi, K0, K1, K2, Vi, Vip1, fip1, fi):
    if K2 == 0.0:
        if K1 == 0.0:
            Kw = (Vip1 - Vi)/K0
        else:
            Kw = np.log(fip1/fi)/K1
    elif Kr < 0.0:
        C1 = (fpip1 - np.sqrt(-Kr))*(fpi + np.sqrt(-Kr))
        C2 = (fpip1 + np.sqrt(-Kr))*(fpi - np.sqrt(-Kr))
        Kw = np.log(C1/C2)/np.sqrt(-Kr)
    elif Kr == 0.0:
        Kw = (2.0/fpi) - (2.0/fpip1)
    elif Kr > 0.0:
        D1 = np.arctan2(fpip1/np.sqrt(Kr))
        D2 = np.arctan2(fpi/np.sqrt(Kr))
        Kw = 2*(D1-D2)/np.sqrt(Kr)
    return Kw


def Kr(K0, K2, K1):
    C1 = 4.0*K0*K2
    C2 = K1**2
    return C1 - C2


def HW3_55(V_w):
    g = 32.2  # ft/sec^2
    S_w = 775  # ft^2
    b_w = 96.4  # ft
    h_w = 11  # ft
    W = 200  # lbf
    mu_r = 0.05
    C_Do = 0.01
    C_Do_L = 0.0
    e = 0.8
    RA = (b_w**2)/S_w
    C_l = 0.4
    C_lmax = 1.4
    Ts = 31.0  # lbf
    Tsp = -1.0  # lbf*sec/ft
    Tspp = 0.0  # lbf*sec^2/ft^2
    t_r = 1  # sec
    rho = 2.3769e-3  # slug/ft^3
    V_LO = VLO(C_lmax, W, S_w, rho)
    Cd = CD(C_Do, C_Do_L, C_l, h_w, b_w, e, RA)
    K_0_1 = K0(Ts, W, mu_r)
    K_1_1 = K1(Tsp, W)
    K_2_1 = K2(Tspp, W, rho, S_w, C_l, mu_r, Cd)
    K_r_1 = Kr(K_0_1, K_2_1, K_1_1)
    fi_1 = fi(K_0_1, K_1_1, V_w, K_2_1)
    fip1_1 = fip1(K_0_1, K_1_1, V_LO, K_2_1)
    fpi_1 = fpi(K_1_1, K_2_1, V_w)
    fpip1_1 = fpi(K_1_1, K_2_1, V_LO)
    K_w_1 = KW(K_r_1, fpip1_1, fpi_1, K_0_1, K_1_1, K_2_1, V_w,
               V_LO, fip1_1, fi_1)
    K_t_1 = KT(V_w, K_0_1, K_1_1, K_2_1, fip1_1, fi_1, V_LO, K_w_1)
    s_a_1 = (K_t_1-(V_w*K_w_1))/g
    s_g = s_a_1 + (V_LO-V_w)*t_r
    print("K0:", K_0_1)
    print("K1:", K_1_1)
    print("K2:", K_2_1)
    print("CD:", Cd)
    print("KR:", K_r_1)
    print("f1:", fi_1)
    print("f2:", fip1_1)
    print("fp1:", fpi_1)
    print("fp2:", fpip1_1)
    print("KW:", K_w_1)
    print("KT:", K_t_1)
    print("Total s_g with HW=",V_w,"ft/s:",s_g)


def HW3_66(V_w):
    g = 32.2  # ft/sec^2
    S_w = 5500  # ft^2
    b_w = 196  # ft
    h_w = 16  # ft
    W = 636600  # lbf
    mu_r = 0.04
    mu_rb = 0.4
    C_Do = 0.033
    C_Do_L = 0.0
    e = 0.64
    RA = (b_w**2)/S_w
    C_l = 0.3
    C_lmax = 1.94
    Ts = 199000  # lbf
    Tsp = -177.6  # lbf*sec/ft
    Tspp = 0.1193  # lbf*sec^2/ft^2
    h_oc = 35  # ft
    t_efr = 1.0  # sec
    t_reac = 2.0  # sec
    t_r = 3.0  # sec
    rho = 2.3769e-3  # slug/ft^3
    Vef = np.arange(1.0, 300.0, 1.0)
    V_LO = VLO(C_lmax, W, S_w, rho)
    V_OC = VOC(C_lmax, W, S_w, rho)
    SOC = np.zeros(len(Vef))
    s_abort = np.zeros(len(Vef))
    Vc = np.zeros(len(Vef))
    for i in range(len(Vef)):
        Cd = CD(C_Do, C_Do_L, C_l, h_w, b_w, e, RA)
        K_0_1 = K0(Ts, W, mu_r)
        K_0_2 = K0(Ts*0.75, W, mu_r)  #
        K_1_1 = K1(Tsp, W)
        K_1_2 = K1(0.75*Tsp, W)  #
        K_2_1 = K2(Tspp, W, rho, S_w, C_l, mu_r, Cd)
        K_2_2 = K2(0.75*Tspp, W, rho, S_w, C_l, mu_r, Cd)  #
        K_r_1 = Kr(K_0_1, K_2_1, K_1_1)
        K_r_2 = Kr(K_0_2, K_2_2, K_1_2)
        fi_1 = fi(K_0_1, K_1_1, V_w, K_2_1)
        fi_2 = fi(K_0_2, K_1_2, Vef[i], K_2_2)
        fip1_1 = fip1(K_0_1, K_1_1, Vef[i], K_2_1)
        fip1_2 = fip1(K_0_2, K_1_2, V_LO, K_2_2)
        fpi_1 = fpi(K_1_1, K_2_1, V_w)
        fpi_2 = fpi(K_1_2, K_2_2, Vef[i])
        fpip1_1 = fpi(K_1_1, K_2_1, Vef[i])
        fpip1_2 = fpi(K_1_2, K_2_2, V_LO)
        K_w_1 = KW(K_r_1, fpip1_1, fpi_1, K_0_1, K_1_1, K_2_1, V_w, Vef[i],
                   fip1_1, fi_1)
        K_w_2 = KW(K_r_2, fpip1_2, fpi_2, K_0_2, K_1_2, K_2_2, Vef[i], V_LO,
                   fip1_2, fi_2)
        K_t_1 = KT(V_w, K_0_1, K_1_1, K_2_1, fip1_1, fi_1, Vef[i], K_w_1)
        K_t_2 = KT(Vef[i], K_0_2, K_1_2, K_2_2, fip1_2, fi_2, V_LO, K_w_2)
        s_a_1 = K_t_1/g
        s_a_2 = K_t_2/g
        s_g = s_a_1 + s_a_2 + V_LO*t_r
        CLOC = (W/(0.5*rho*(V_OC**2)*S_w))
        CDOC = C_Do + (C_Do_L*CLOC) + ((CLOC**2)/(np.pi*e*RA))
        DOC = 0.5*rho*(V_OC**2)*S_w*CDOC
        G_OC = np.arcsin(((0.75*Ts) - DOC)/W)  #
        CLOC2 = W*np.cos(G_OC)/(0.5*rho*(V_OC**2)*S_w)
        CDOC2 = CD(C_Do, C_Do_L, CLOC2, (h_w+h_oc), b_w, e, RA)
        DOC2 = 0.5*rho*(V_OC**2)*S_w*CDOC2
        G_OC2 = np.arcsin(((0.75*T(V_OC)) - DOC2)/W)
        CLLO = W/(0.5*rho*(V_LO**2)*S_w)
        CDLO = CD(C_Do, C_Do_L, CLLO, h_w, b_w, e, RA)
        DLO = 0.5*rho*(V_LO**2)*S_w*CDLO
        Fc = 0.5*(0.75*T(V_LO) - DLO +((0.75*T(V_OC) - DOC2)/np.cos(G_OC2)))  #
        s_c_C1 = ((V_OC**2)-V_LO**2)/(2*g)
        s_c = (W/Fc)*(h_oc + s_c_C1)
        SOC[i] = s_g + s_c
        V_fr = Vef[i] + g*(K_0_2 + K_1_2*Vef[i] + K_2_2*(Vef[i]**2))*(t_reac - t_efr)
        V_frr = V_fr + g*(K_0_2 + K_1_2*V_fr + K_2_2*(V_fr**2))*(t_reac - t_efr)
        fi_r = fi(K_0_2, K_1_2, Vef[i], K_2_2)
        fipr = fip1(K_0_2, K_1_2, V_frr, K_2_2)
        fpir = fpi(K_1_2, K_2_2, Vef[i])
        fpipr = fpi(K_1_2, K_2_2, V_frr)
        K_w_r = KW(K_r_2, fpipr, fpir, K_0_2, K_1_2, K_2_2, Vef[i], V_frr,
                   fpir, fi_r)
        K_t_r = KT(Vef[i], K_0_2, K_1_2, K_2_2, fipr, fi_r, V_frr, K_w_r)
        s_frr = K_t_r/g
        K_0_b = K0(V_w, W, mu_rb)
        K_1_b = K1(V_w, W)
        K_2_b = K2(V_w, W, rho, S_w, C_l, mu_rb, Cd)
        K_r_b = Kr(K_0_b, K_2_b, K_1_b)
        fi_b = fi(K_0_b, K_1_b, V_frr, K_2_b)
        fipb = fip1(K_0_b, K_1_b, V_w, K_2_b)
        fpib = fpi(K_1_b, K_2_b, V_frr)
        fpipb = fpi(K_1_b, K_2_b, V_w)
        K_w_b = KW(K_r_b, fpipb, fpib, K_0_b, K_1_b, K_2_b, V_frr, V_w,
                   fpib, fi_b)
        K_t_b = KT(V_frr, K_0_b, K_1_b, K_2_b, fipb, fi_b, V_w, K_w_b)
        sb = (K_t_b/g)
#        sb2 = ((W/S_w)/(rho*g*(Cd - mu_r*C_l)))*np.log(1. + (rho*(V_frr*V_frr)/(2.*W/S_w))*(Cd/mu_r - C_l))
        s_abort[i] = s_a_1 + s_frr + sb
        if Vef[i] == 50.0:
            print("s_OC at V_ef = 50ft/s: ",SOC[i])
        plt.figure(1)
        plt.scatter(SOC, Vef)
        plt.scatter(s_abort, Vef)
    for i in range(len(Vef)):
        Vc[i] = np.abs(SOC[i] - s_abort[i])
    Vcefr = Vef[np.argmin(Vc)]
    bfl = SOC[np.argmin(Vc)]
    print("Critical engine failure speed : ",Vcefr,
          " ft/s\nBalanced Field Length = ", bfl, "ft")

print("3.55")
print("-----------------------")
HW3_55(0.*1.467)
HW3_55(5.*1.467)
print("\n3.66")
print("-----------------------")
HW3_66(0.0*1.467)
