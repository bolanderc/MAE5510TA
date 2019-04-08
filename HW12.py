#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 09:09:26 2019

@author: christian
"""

import numpy as np
import scipy.linalg
from eig_vec_norm import normalize_evec
import matplotlib.pyplot as plt
from matplotlib import rc



rc('text', usetex=True)

def HW9_46():
    Sw = 185
    bw = 33
    cw = Sw/bw
    W = 2800
    rho = 0.0023769
    g = 32.17
    Vo = 180
    M = 0.0
    Ixxb = 1000
    Iyyb = 3000
    Izzb = 3500
    Ixyb = 0.0
    Ixzb = 30
    Iyzb = 0
    hxb = 0.
    hyb = 0.
    hzb = 0.
    CLo = 0.393
    CDo = 0.05
    Cmo = 0.
    CLuhat = 0.
    CDuhat = 0.
    Cmuhat = 0.
    CLahat = 1.6
    CDahat = 0.
    Cmahat = -4.35
    CLM = 0.
    CDM = 0.
    CmM = 0.
    CLa = 4.4
    CDa = 0.35
    Cella = 0.
    Cma = -0.68
    Cna = 0.
    CYb = -0.56
    Cellb = -0.075
    Cmb = 0.
    Cnb = 0.07
    CYp = 0.
    Cellp = -0.41
    Cnp = -0.0575
    CLq = 3.8
    CDq = 0.
    Cmq = -9.95
    CYr = 0.24
    Cellr = 0.105
    Cnr = -0.125
    theta_o = 0.
    phi_o = 60.
    Cto = np.cos(theta_o)
    Cpo = np.cos(phi_o*np.pi/180)
    CLref = W/(0.5*rho*Vo*Vo*Sw)
    CL = W*Cto/(0.5*rho*Vo*Vo*Sw*Cpo)
    CD = CDo + (CDa*CLref/(2.*CLa))*((CL*CL/(CLref*CLref)) - 1.)
    CDa = CDa*CL/CLref
    CXa = CL - CDa
    CZa = -CLa - CD
    delta_a = (CL - CLref)/CLa
    phi = -delta_a
    Cp = np.cos(phi)
    Sp = np.sin(phi)
    Ixx_s = Ixxb*Cp*Cp + 2.*Ixzb*Cp*Sp + Izzb*Sp*Sp
    Iyy_s = Iyyb
    Izz_s = Izzb*Cp*Cp - 2.*Ixzb*Cp*Sp + Ixxb*Sp*Sp
    Ixz_s = Ixzb*(Cp*Cp - Sp*Sp) + (Izzb - Ixxb)*Cp*Sp
    Ixy_s = Ixyb
    Iyz_s = Iyzb
    CX = -CD - M*CDM
    CZ = -CL - M*CLM
    Cm = Cmo
    CXuhat = -CDuhat
    CZuhat = -CLuhat
    CXahat = -CDahat
    CZahat = -CLahat
    CXuhat_s = CXuhat*Cp*Cp - (CXahat + CZuhat)*Cp*Sp + CZahat*Sp*Sp
    CZuhat_s = CXahat*Cp*Cp + (CXuhat - CZahat)*Cp*Sp - CZuhat*Sp*Sp
    Cmuhat_s = Cmuhat*Cp - Cmahat*Sp
    CXahat_s = CXahat*Cp*Cp + (CXuhat - CZahat)*Cp*Sp - CZuhat*Sp*Sp
    CZahat_s = CZahat*Cp*Cp + (CXahat + CZuhat)*Cp*Sp + CXuhat*Sp*Sp
    Cmahat_s = Cmahat*Cp + Cmuhat*Sp
    CXu_s = -M*CDM
    CZu_s = -M*CLM
    Cmu_s = -M*CmM
    CXa_s = CXa
    CZa_s = CZa
    Cma_s = Cma*Cp + Cmu_s*Sp
    CYb_s = CYb
    Cellb_s = Cellb*Cp - Cnb*Sp
    Cnb_s = Cnb*Cp + Cellb*Sp
    Cella_s = Cella
    Cmb_s = Cmb
    Cna_s = Cna
    CYp_s = CYp*Cp - CYr*Sp
    Cellp_s = Cellp*Cp*Cp - (Cellr + Cnp)*Cp*Sp + Cnr*Sp*Sp
    Cnp_s = Cnp*Cp*Cp + (Cellp - Cnr)*Cp*Sp - Cellr*Sp*Sp
    CXq = -CDq
    CZq = -CLq
    CXq_s = CXq*Cp - CZq*Sp
    CZq_s = CZq*Cp + CXq*Sp
    Cmq_s = Cmq
    CYr_s = CYr*Cp + CYp*Sp
    Cellr_s = Cellr*Cp*Cp + (Cellp - Cnr)*Cp*Sp - Cnp*Sp*Sp
    Cnr_s = Cnr*Cp*Cp + (Cellr + Cnp)*Cp*Sp + Cellp*Sp*Sp
    l_ref = 13.6
    ixy = Ixy_s/Ixx_s
    ixz = Ixz_s/Ixx_s
    iyx = Ixy_s/Iyy_s
    iyz = Iyz_s/Iyy_s
    izx = Ixz_s/Izz_s
    izy = Iyz_s/Izz_s
    Bm = rho*Sw*cw/(4.*W/g)
    By = rho*Sw*cw*cw*l_ref/(4.*Iyy_s)
    Spo = np.sin(phi_o*np.pi/180)
    Tpo = Spo/Cpo
    Sto = np.sin(theta_o*np.pi/180)
    Tto = Sto/Cto
    Bxup = Bm*CXuhat_s
    Bzup = Bm*CZuhat_s
    Bmup = By*Cmuhat_s
    Bxap = Bm*CXahat_s
    Bzap = Bm*CZahat_s
    Bmap = By*Cmahat_s
    Ag = g*l_ref/(Vo*Vo)
    Am = rho*Sw*l_ref/(2.*W/g)
    Ay = rho*Sw*cw*l_ref*l_ref/(2.*Iyy_s)
    etaxx = Ag*(Ixz_s*Tpo*Spo*Cto - Ixy_s*Spo*Cto)/Ixx_s
    etaxy = hzb*l_ref/(Ixx_s*Vo) + Ag*((Izz_s - Iyy_s)*Spo*Cto -
                                       2.*Iyz_s*Tpo*Spo*Cto +
                                       Ixz_s*Tpo*Sto)/Ixx_s
    etaxz = hyb*l_ref/(Ixx_s*Vo) + Ag*((Iyy_s - Izz_s)*Tpo*Spo*Cto -
                                       2.*Iyz_s*Spo*Cto +
                                       Ixy_s*Tpo*Sto)/Ixx_s
    etayx = hzb*l_ref/(Iyy_s*Vo) + Ag*((Izz_s - Ixx_s)*Spo*Cto +
                                       2.*Ixz_s*Tpo*Sto -
                                       Iyz_s*Tpo*Spo*Cto)/Iyy_s
    etayy = Ag*(Ixy_s*Spo*Cto + Iyz_s*Tpo*Sto)/Iyy_s
    etayz = hxb*l_ref/(Iyy_s*Vo) + Ag*((Izz_s - Ixx_s)*Tpo*Sto -
                                       2.*Ixz_s*Spo*Cto -
                                       Ixy_s*Tpo*Spo*Cto)/Iyy_s
    etazx = hyb*l_ref/(Izz_s*Vo) + Ag*((Iyy_s - Ixx_s)*Tpo*Spo*Cto +
                                       2.*Ixy_s*Tpo*Sto - Iyz_s*Spo*Cto)/Izz_s
    etazy = hxb*l_ref/(Izz_s*Vo) + Ag*((Iyy_s - Ixx_s)*Tpo*Sto -
                                       2.*Ixy_s*Tpo*Spo*Cto -
                                       Ixz_s*Spo*Cto)/Izz_s
    etazz = Ag*(-Iyz_s*Tpo*Sto - Ixz_s*Tpo*Spo*Cto)/Izz_s
    Axu = Am*(2.*CX + CXu_s)
    Azu = Am*(2.*CZ + CZu_s)
    Amu = Ay*(2.*Cm + Cmu_s)
    Amlong = rho*Sw*cw/(4.*W/g)
    Amlat = rho*Sw*bw/(4.*W/g)
    Axa = Am*CXa_s
    Aza = Am*CZa_s
    Ama = Ay*Cma_s
    Ayb = Am*CYb_s
    Ax = rho*Sw*bw*l_ref*l_ref/(2.*Ixx_s)
    Aellb = Ax*Cellb_s
    Az = rho*Sw*bw*l_ref*l_ref/(2.*Izz_s)
    Anb = Az*Cnb_s
    Aella = Ax*Cella_s
    Amb = Ay*Cmb_s
    Ana = Az*Cna_s
    Ayp = Amlat*CYp_s
    Axlat = rho*Sw*bw*bw*l_ref/(4.*Ixx_s)
    Aellp = Axlat*Cellp_s
    Azlat = rho*Sw*bw*bw*l_ref/(4.*Izz_s)
    Anp = Azlat*Cnp_s
    Axq = Amlong*CXq_s
    Azq = Amlong*CZq_s
    Aylong = rho*Sw*cw*cw*l_ref/(4.*Iyy_s)
    Amq = Aylong*Cmq_s
    Ayr = Amlat*CYr_s
    Aellr = Axlat*Cellr_s
    Anr = Azlat*Cnr_s
    A = np.array([[Axu, Ag*Spo*Cto, Axa - Ag*Tpo*Spo*Cto, 0., Axq, 0., 0., 0., 0., 0., -Ag*Cto, 0.],
                  [-Ag*Spo*Cto, Ayb, -Ag*Tpo*Sto, Ayp, 0., Ayr - 1., 0., 0., 0., Ag*Cpo*Cto, -Ag*Spo*Sto, 0.],
                  [Azu + Ag*Tpo*Spo*Cto, Ag*Tpo*Sto, Aza, 0., Azq + 1., 0., 0., 0., 0., -Ag*Spo*Cto, -Ag*Cpo*Sto, 0.],
                  [0., Aellb, Aella, Aellp + etaxx, -etaxy, Aellr + etaxz, 0., 0., 0., 0., 0., 0.],
                  [Amu, Amb, Ama, etayx, Amq + etayy, -etayz, 0., 0., 0., 0., 0., 0.],
                  [0., Anb, Ana, Anp - etazx, etazy, Anr + etazz, 0., 0., 0., 0., 0., 0.],
                  [Cto, Spo*Sto, Cpo*Sto, 0., 0., 0., 0., 0., 0., 0., -Sto, 0.],
                  [0., Cpo, -Spo, 0., 0., 0., 0., 0., 0., 0., 0., Cto],
                  [-Sto, Spo*Cto, Cpo*Cto, 0., 0., 0., 0., 0., 0., 0., -Cto, 0.],
                  [0., 0., 0., 1., Spo*Tto, Cpo*Tto, 0., 0., 0., 0., Ag*Tpo/Cto, 0.],
                  [0., 0., 0., 0., Cpo, -Spo, 0., 0., 0., -Ag*Tpo*Cto, 0., 0.],
                  [0., 0., 0., 0., Spo/Cto, Cpo/Cto, 0., 0., 0., 0., Ag*Tpo*Tto, 0.]])
    B = np.zeros((12, 12))
    np.fill_diagonal(B, 1.)
    B[0, 0] = 1. - Bxup
    B[0, 2] = -Bxap
    B[2, 0] = -Bzup
    B[2, 2] = 1 - Bzap
    B[3, 4:6] = np.array([-ixy, -ixz])
    B[4, 0:4] = np.array([-Bmup, 0., -Bmap, iyx])
    B[4, 5] = -iyz
    B[5, 3:5] = np.array([-izx, -izy])
    eigvals, eigvecs = scipy.linalg.eig(A, B)
    eigv_norm = normalize_evec(eigvecs)
    rigid = []
    dutch = []
    short = []
    long = []
    osc = []
    maxr = np.argmax(abs(eigvals.real))
    for i in range(len(eigvals)):
        if eigvals[i].imag != 0.0:
            osc.append(i)
    maxi = osc[np.argmax(abs(eigvals[osc].real))]
    mini = osc[np.argmin(abs(eigvals[osc].real))]
    for i in range(len(eigvals)):
        if eigvals[i].real == 0.0 and eigvals[i].imag == 0.0:
            rigid.append(i)
        elif eigvals[i].imag == 0.0:
            if np.round(eigvals[i].real, decimals=6) == np.round(eigvals[maxr].real, decimals=6):
                roll = i
            else:
                spiral = i
        else:
            if np.round(eigvals[i].real, decimals=6) == np.round(eigvals[maxi].real, decimals=6):
                short.append(i)
            elif np.round(eigvals[i].real, decimals=6) == np.round(eigvals[mini].real, decimals=6):
                long.append(i)
            else:
                dutch.append(i)
    names = [r'$\Delta\mu$',  r'$\Delta\beta$', r'$\Delta\alpha$', r"$\Delta p$",
             r"$\Delta q$", r"$\Delta r$", r'$\Delta\varsigma_x$', r'$\Delta\varsigma_y$',
             r'$\Delta\varsigma_z$', r'$\Delta\phi$',  r'$\Delta\theta$',
             r'$\Delta\psi$']
    for i in range(len(eigvals)):
        if i in rigid:
            plotname='Rigid-body'
            print("Rigid-body Mode")
            print("----------------")
            print(u"\u03BB" + ' =  (%7.6f + %7.6fj)' % (eigvals[i].real,
                                                        eigvals[i].imag))
        elif i == roll:
            plotname = 'Roll'
            print("Roll Mode")
            print("----------------")
            print(u"\u03BB" + ' =  (%7.6f + %7.6fj)' % (eigvals[i].real,
                                                        eigvals[i].imag))
            sqrtmult = np.sqrt(eigvals[roll]*eigvals[spiral])
            sigma = -eigvals[i].real*2.*Vo/bw
            dampingtime = np.log(0.01)/-sigma
            wn = 2.*Vo*sqrtmult/bw
            damp_ratio = -1.*(eigvals[roll] + eigvals[spiral])/(2.*sqrtmult)
            print(u"\u03C3" + " = %7.6f sec^-1" % (sigma))
            print("99 %% damping time = %3.2f sec" % (dampingtime))
            print(u"\u03C9_n" + " = %7.6f rad/sec" % wn)
            print(u"\u03B6" + " = %7.6f" % damp_ratio)
        elif i in dutch:
            plotname = 'Dutch Roll'
            print("Dutch Roll Mode")
            print("----------------")
            print(u"\u03BB" + ' =  (%7.6f + %7.6fj)' % (eigvals[i].real,
                                                        eigvals[i].imag))
            sqrtmult = np.sqrt(eigvals[dutch[0]]*eigvals[dutch[1]])
            sigma = -eigvals[i].real*2.*Vo/bw
            dampingtime = np.log(0.01)/-sigma
            wd = 2.*Vo*np.abs(eigvals[i].imag)/bw
            period = 2.*np.pi/wd
            wn = 2.*Vo*sqrtmult/bw
            damp_ratio = -1.*(eigvals[dutch[0]] + eigvals[dutch[1]])/(2.*sqrtmult)
            print(u"\u03C3" + " = %7.6f sec^-1" % (sigma))
            print("99 %% damping time = %3.2f sec" % (dampingtime))
            print(u"\u03C9_d" + " = %7.6f rad/sec" % wd)
            print("period = %7.6f" % period)
            print(u"\u03C9_n" + " = %7.6f rad/sec" % wn)
            print(u"\u03B6" + " = %7.6f" % damp_ratio)
        elif i == spiral:
            plotname = 'Spiral'
            print("Spiral Mode")
            print("----------------")
            print(u"\u03BB" + ' =  (%7.6f + %7.6fj)' % (eigvals[i].real,
                                                        eigvals[i].imag))
            sqrtmult = np.sqrt(eigvals[roll]*eigvals[spiral])
            sigma = -eigvals[i].real*2.*Vo/bw
            dampingtime = np.log(0.01)/-sigma
            wn = 2.*Vo*sqrtmult/bw
            damp_ratio = -1.*(eigvals[roll] + eigvals[spiral])/(2.*sqrtmult)
            print(u"\u03C3" + " = %7.6f sec^-1" % (sigma))
            print("99 %% damping time = %3.2f sec" % (dampingtime))
            print(u"\u03C9_n" + " = %7.6f rad/sec" % wn)
            print(u"\u03B6" + " = %7.6f" % damp_ratio)
        elif i in long:
            plotname = 'Long Period (Phugoid)'
            print("Long Period (Phugoid) Mode")
            print("----------------")
            print(u"\u03BB" + ' =  (%7.6f + %7.6fj)' % (eigvals[i].real,
                                                        eigvals[i].imag))
            sqrtmult = np.sqrt(eigvals[long[0]]*eigvals[long[1]])
            sigma = -eigvals[i].real*2.*Vo/bw
            dampingtime = np.log(0.01)/-sigma
            wd = 2.*Vo*np.abs(eigvals[i].imag)/bw
            period = 2.*np.pi/wd
            wn = 2.*Vo*sqrtmult/bw
            damp_ratio = -1.*(eigvals[long[0]] + eigvals[long[1]])/(2.*sqrtmult)
            print(u"\u03C3" + " = %7.6f sec^-1" % (sigma))
            print("99 %% damping time = %3.2f sec" % (dampingtime))
            print(u"\u03C9_d" + " = %7.6f rad/sec" % wd)
            print("period = %7.6f" % period)
            print(u"\u03C9_n" + " = %7.6f rad/sec" % wn)
            print(u"\u03B6" + " = %7.6f" % damp_ratio)
        elif i in short:
            plotname = 'Short Period'
            print("Short Period Mode")
            print("----------------")
            print(u"\u03BB" + ' =  (%7.6f + %7.6fj)' % (eigvals[i].real,
                                                        eigvals[i].imag))
            sqrtmult = np.sqrt(eigvals[short[0]]*eigvals[short[1]])
            sigma = -eigvals[i].real*2.*Vo/bw
            dampingtime = np.log(0.01)/-sigma
            wd = 2.*Vo*np.abs(eigvals[i].imag)/bw
            period = 2.*np.pi/wd
            wn = 2.*Vo*sqrtmult/bw
            damp_ratio = -1.*(eigvals[short[0]] + eigvals[short[1]])/(2.*sqrtmult)
            print(u"\u03C3" + " = %7.6f sec^-1" % (sigma))
            print("99 %% damping time = %3.2f sec" % (dampingtime))
            print(u"\u03C9_d" + " = %7.6f rad/sec" % wd)
            print("period = %7.6f" % period)
            print(u"\u03C9_n" + " = %7.6f rad/sec" % wn)
            print(u"\u03B6" + " = %7.6f" % damp_ratio)
        print(u"\u03C7" + ' =  (%7.6f + %7.6fj)' % (eigv_norm[i, 0].real,
                                                    eigv_norm[i, 0].imag))
        t = np.arange(0, 100, 0.01)
        y = eigv_norm[i, 0]*np.exp(eigvals[i]*t)
        plt.figure(i, figsize = (6, 6))
        plt.plot(t, y, label=names[0])
        for j in range(len(eigvals) - 1):
            print('     (%7.6f + %7.6fj)' % (eigv_norm[i, j + 1].real,
                                             eigv_norm[i, j + 1].imag))
            y = eigv_norm[i, j + 1]*np.exp(eigvals[i]*t)
            plt.plot(t, y, label=names[j + 1])
        plt.legend(loc='upper right')
        plt.title(plotname)
        plt.show()







HW9_46()
