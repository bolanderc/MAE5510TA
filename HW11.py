#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 10:09:58 2019

@author: christian
"""

import numpy as np
import scipy.linalg
from eig_vec_norm import normalize_evec

import warnings
warnings.filterwarnings('ignore')


def matprint(mat, fmt="g"):
    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
        print("")

def HW8_3():
    rho = 0.0023769
    g = 32.17
    theta_o = 0.0
    Sw = 185
    bw = 33
    W = 2800
    Vo = 180
    CDo = 0.05
    Ixxb = 1000
    Iyyb = 3000
    Izzb = 3500
    Ixzb = 30
    CLa = 4.4
    CDa = 0.35
    Cma = -0.68
    CLahat = 1.6
    Cmahat = -4.35
    CYb = -0.56
    Clb = -0.075
    Cnb = 0.07
    CDq = 0.
    CLq = 3.8
    Cmq = -9.95
    CYp = 0.
    Clp = -0.410
    Cnp = -0.0575
    CYr = 0.24
    Clr = 0.105
    Cnr = -0.125
    W_g = W/g
    cw = Sw/bw
    Cg = rho*Sw*cw/(4*W_g)
    Cy = rho*Sw*cw*cw*cw/(8.*Iyyb)
    Rgx = g*cw/(2.*Vo*Vo)
    CLo = W*np.cos(theta_o)/(0.5*rho*Vo*Vo*Sw)
    Rzahat = Cg*-CLahat
    Rmahat = Cy*Cmahat
    Rxmu = Cg*2.*-CDo
    Rzmu = -Cg*2.*CLo
    Rmmu = 0.0
    Rxa = Cg*(CLo - CDa)
    Rza = Cg*(-CLa - CDo)
    Rma = Cy*Cma
    Rxq = -Cg*CDq
    Rzq = -Cg*CLq
    Rmq = Cy*Cmq
    Ct = np.cos(theta_o)
    St = np.sin(theta_o)
    A = np.array([[Rxmu, Rxa, Rxq, 0., 0., -Rgx*Ct],
                  [Rzmu, Rza, (1 + Rzq), 0., 0., -Rgx*St],
                  [Rmmu, Rma, Rmq, 0., 0., 0.],
                  [Ct, St, 0., 0., 0., -St],
                  [-St, Ct, 0., 0., 0., -Ct],
                  [0., 0., 1., 0., 0., 0.]])
    B = np.array([[1., 0., 0., 0., 0., 0.],
                  [0., (1. - Rzahat), 0., 0., 0., 0.],
                  [0., -Rmahat, 1., 0., 0., 0.],
                  [0., 0., 0., 1., 0., 0.],
                  [0., 0., 0., 0., 1., 0.],
                  [0., 0., 0., 0., 0., 1.]])
    eigvals, eigvecs = scipy.linalg.eig(A, B)
    eigv_norm = normalize_evec(eigvecs)
    for i in range(len(eigvecs)):
        print('eig_' + str(i + 1) + '    =  (%7.6f + %7.6fj)' % (eigvals[i].real,
                                                                 eigvals[i].imag))
        print('vec      =  (%7.6f + %7.6fj)' % (eigvecs[0, i].real,
                                                eigvecs[0, i].imag))
        print('            (%7.6f + %7.6fj)' % (eigvecs[1, i].real,
                                                eigvecs[1, i].imag))
        print('            (%7.6f + %7.6fj)' % (eigvecs[2, i].real,
                                                eigvecs[2, i].imag))
        print('            (%7.6f + %7.6fj)' % (eigvecs[3, i].real,
                                                eigvecs[3, i].imag))
        print('            (%7.6f + %7.6fj)' % (eigvecs[4, i].real,
                                                eigvecs[4, i].imag))
        print('            (%7.6f + %7.6fj)' % (eigvecs[5, i].real,
                                                eigvecs[5, i].imag))
        print('vec_norm =  (%7.6f + %7.6fj)' % (eigv_norm[i, 0].real,
                                                eigv_norm[i, 0].imag))
        print('            (%7.6f + %7.6fj)' % (eigv_norm[i, 1].real,
                                                eigv_norm[i, 1].imag))
        print('            (%7.6f + %7.6fj)' % (eigv_norm[i, 2].real,
                                                eigv_norm[i, 2].imag))
        print('            (%7.6f + %7.6fj)' % (eigv_norm[i, 3].real,
                                                  eigv_norm[i, 3].imag))
        print('            (%7.6f + %7.6fj)' % (eigv_norm[i, 4].real,
                                                  eigv_norm[i, 4].imag))
        print('            (%7.6f + %7.6fj)\n' % (eigv_norm[i, 5].real,
                                                  eigv_norm[i, 5].imag))


def HW8_4():
    Ct = np.cos(0.)
    St = np.sin(0.)
    cw = 28.
    Vo = 660.
    Rgx = 0.00103
    Rxmuhat = 0.
    Rzmuhat = 0.
    Rmmuhat = 0.
    Rxahat = 0.
    Rzahat = -0.00239
    Rmahat = -0.00201
    Rxmu = -0.000206
    Rzmu = -0.00288
    Rmmu = 0.
    Rxa = 0.000354
    Rza = -0.00960
    Rma = -0.000512
    Rxq = 0.
    Rzq = -0.0102
    Rmq = -0.00857
    A = np.array([[Rxmu, Rxa, Rxq, 0., 0., -Rgx*Ct],
                  [Rzmu, Rza, (1 + Rzq), 0., 0., -Rgx*St],
                  [Rmmu, Rma, Rmq, 0., 0., 0.],
                  [Ct, St, 0., 0., 0., -St],
                  [-St, Ct, 0., 0., 0., -Ct],
                  [0., 0., 1., 0., 0., 0.]])
    B = np.array([[1., 0., 0., 0., 0., 0.],
                  [0., (1. - Rzahat), 0., 0., 0., 0.],
                  [0., -Rmahat, 1., 0., 0., 0.],
                  [0., 0., 0., 1., 0., 0.],
                  [0., 0., 0., 0., 1., 0.],
                  [0., 0., 0., 0., 0., 1.]])
    eigvals, eigvecs = scipy.linalg.eig(A, B)
    eigv_norm = normalize_evec(eigvecs)
    sp = np.argmax(np.absolute(eigvals))
    for i in range(len(eigvals)):
        if eigvals[i].real == 0.0 and eigvals[i].imag == 0.0:
            print("Rigid-body Mode")
            print("----------------")
            print(u"\u03BB" + ' =  (%7.6f + %7.6fj)' % (eigvals[i].real,
                                                        eigvals[i].imag))
        elif np.round(eigvals[i].real, decimals=6) == np.round(eigvals[sp].real, decimals=6):
            print("Short-period Mode")
            print("----------------")
            print(u"\u03BB" + ' =  (%7.6f + %7.6fj)' % (eigvals[i].real,
                                                        eigvals[i].imag))
            dampingtime = np.log(0.01)*cw/(eigvals[i].real*2.*Vo)
            period = 2.*np.pi*cw/(np.abs(eigvals[i].imag)*2.*Vo)
            print("99 %% damping time = %3.2f sec" % (dampingtime))
            print("period = %3.2f sec" % (period))
        else:
            print("Long-period (Phugoid) Mode")
            print("----------------")
            print(u"\u03BB" + ' =  (%7.6f + %7.6fj)' % (eigvals[i].real,
                                                        eigvals[i].imag))
            dampingtime = np.log(0.01)*cw/(eigvals[i].real*2.*Vo)
            period = 2.*np.pi*cw/(np.abs(eigvals[i].imag)*2.*Vo)
            print("99 %% damping time = %3.2f sec" % (dampingtime))
            print("period = %3.2f sec" % (period))
        print(u"\u03C7" + ' =  (%7.6f + %7.6fj)' % (eigv_norm[i, 0].real,
                                                    eigv_norm[i, 0].imag))
        print('     (%7.6f + %7.6fj)' % (eigv_norm[i, 1].real,
                                         eigv_norm[i, 1].imag))
        print('     (%7.6f + %7.6fj)' % (eigv_norm[i, 2].real,
                                         eigv_norm[i, 2].imag))
        print('     (%7.6f + %7.6fj)' % (eigv_norm[i, 3].real,
                                         eigv_norm[i, 3].imag))
        print('     (%7.6f + %7.6fj)' % (eigv_norm[i, 4].real,
                                         eigv_norm[i, 4].imag))
        print('     (%7.6f + %7.6fj)\n' % (eigv_norm[i, 5].real,
                                           eigv_norm[i, 5].imag))


def HW9_1():
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
    A = np.array([[RyB, Ryp, (Ryr - 1), 0., Rgy*np.cos(theta_o), 0.,],
                  [RlB, Rlp, Rlr, 0., 0., 0.,],
                  [RnB, Rnp, Rnr, 0., 0., 0.,],
                  [1., 0., 0., 0., 0., np.cos(theta_o)],
                  [0., 1., np.tan(theta_o), 0., 0., 0.],
                  [0., 0., (1./np.cos(theta_o)), 0., 0., 0.]])
    B = np.array([[1., 0., 0., 0., 0., 0.],
                  [0., 1., -ixz, 0., 0., 0.],
                  [0., -izx, 1., 0., 0., 0.],
                  [0., 0., 0., 1., 0., 0.],
                  [0., 0., 0., 0., 1., 0.],
                  [0., 0., 0., 0., 0., 1.]])
    eigvals, eigvecs = scipy.linalg.eig(A, B)
    matprint(B)
    eigv_norm = normalize_evec(eigvecs)
    sp = np.argmax(abs(eigvals.real))
    rigid = []
    dutch = []
    for i in range(len(eigvals)):
        if eigvals[i].real == 0.0 and eigvals[i].imag == 0.0:
            rigid.append(i)
        elif eigvals[i].imag != 0.0:
            dutch.append(i)
        elif np.round(eigvals[i].real, decimals=6) == np.round(eigvals[sp].real, decimals=6):
            roll = i
        else:
            spiral = i
    for i in range(len(eigvals)):
        if i in rigid:
            print("Rigid-body Mode")
            print("----------------")
            print(u"\u03BB" + ' =  (%7.6f + %7.6fj)' % (eigvals[i].real,
                                                        eigvals[i].imag))
        elif i == roll:
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
        else:
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
        print(u"\u03C7" + ' =  (%7.6f + %7.6fj)' % (eigv_norm[i, 0].real,
                                                    eigv_norm[i, 0].imag))
        print('     (%7.6f + %7.6fj)' % (eigv_norm[i, 1].real,
                                         eigv_norm[i, 1].imag))
        print('     (%7.6f + %7.6fj)' % (eigv_norm[i, 2].real,
                                         eigv_norm[i, 2].imag))
        print('     (%7.6f + %7.6fj)' % (eigv_norm[i, 3].real,
                                         eigv_norm[i, 3].imag))
        print('     (%7.6f + %7.6fj)' % (eigv_norm[i, 4].real,
                                         eigv_norm[i, 4].imag))
        print('     (%7.6f + %7.6fj)\n' % (eigv_norm[i, 5].real,
                                           eigv_norm[i, 5].imag))
    rollamp = np.absolute(eigv_norm[roll, :])
    spiralamp = np.absolute(eigv_norm[spiral, :])
    dutchamp = np.absolute(eigv_norm[dutch[0], :])
    dutchphase = np.arctan2(eigv_norm[dutch[0], :].imag, eigv_norm[dutch[0], :].real)*180/np.pi
    print('Eigenvector\tRoll Mode\tSpiral Mode\t\tDutch Roll')
    print('Component\tAmplitude\tAmplitude\tAmplitude\tPhase')
    print('----------------------------------------------------------------------')
    symbols = [u"\u0392", "p", "r", u"\u03BEy", u"\u03D5", u"\u03C8"]
    for i in range(len(rollamp)):
        print('    ' + u"\u0394" + '%s\t\t%7.6f\t%7.6f\t%7.6f     %4.2f\u00B0' % (symbols[i],
                                                                                    rollamp[i],
                                                                                    spiralamp[i],
                                                                                    dutchamp[i],
                                                                                    dutchphase[i]))


def HWEx9_2_1():
    Sw = 185
    bw = 33
    W = 2800
    Vo = 180
    CDo = 0.05
    Ixxb = 1000
    Iyyb = 3000
    Izzb = 3500
    Ixzb = 30
    CLa = 4.4
    CDa = 0.35
    Cma = -0.68
    CLahat = 1.6
    Cmahat = -4.35
    CYb = -0.56
    Clb = -0.075
    Cnb = 0.07
    CDq = 0.
    CLq = 3.8
    Cmq = -9.95
    CYp = 0.
    Clp = -0.41
    Cnp = -0.0575
    CYr = 0.24
    Clr = 0.105
    Cnr = -0.125
    ixz = Ixzb/Ixxb
    izx = Ixzb/Izzb
    Cm = 0.041679
    Cx = 1.975306
    Cz = 0.564373
    RyB = Cm*CYb
    RlB = Cx*Clb
    RnB = Cz*Cnb
    Ryp = Cm*CYp
    Rlp = Cx*Clp
    Rnp = Cz*Cnp
    Ryr = Cm*CYr
    Rlr = Cx*Clr
    Rnr = Cz*Cnr
    Rgy = 0.0164
    theta_o = 0.0
    A = np.array([[RyB, Ryp, (Ryr - 1), 0., Rgy*np.cos(theta_o), 0.,],
                  [RlB, Rlp, Rlr, 0., 0., 0.,],
                  [RnB, Rnp, Rnr, 0., 0., 0.,],
                  [1., 0., 0., 0., 0., np.cos(theta_o)],
                  [0., 1., np.tan(theta_o), 0., 0., 0.],
                  [0., 0., (1./np.cos(theta_o)), 0., 0., 0.]])
    B = np.array([[1., 0., 0., 0., 0., 0.],
                  [0., 1., -ixz, 0., 0., 0.],
                  [0., -izx, 1., 0., 0., 0.],
                  [0., 0., 0., 1., 0., 0.],
                  [0., 0., 0., 0., 1., 0.],
                  [0., 0., 0., 0., 0., 1.]])
    eigvals, eigvecs = scipy.linalg.eig(A, B)
    eigv_norm = normalize_evec(eigvecs)
    sp = np.argmax(abs(eigvals.real))
    rigid = []
    dutch = []
    for i in range(len(eigvals)):
        if eigvals[i].real == 0.0 and eigvals[i].imag == 0.0:
            rigid.append(i)
        elif eigvals[i].imag != 0.0:
            dutch.append(i)
        elif np.round(eigvals[i].real, decimals=6) == np.round(eigvals[sp].real, decimals=6):
            roll = i
        else:
            spiral = i
    for i in range(len(eigvals)):
        if i in rigid:
            print("Rigid-body Mode")
            print("----------------")
            print(u"\u03BB" + ' =  (%7.6f + %7.6fj)' % (eigvals[i].real,
                                                        eigvals[i].imag))
        elif i == roll:
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
        else:
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
        print(u"\u03C7" + ' =  (%7.6f + %7.6fj)' % (eigv_norm[i, 0].real,
                                                    eigv_norm[i, 0].imag))
        print('     (%7.6f + %7.6fj)' % (eigv_norm[i, 1].real,
                                         eigv_norm[i, 1].imag))
        print('     (%7.6f + %7.6fj)' % (eigv_norm[i, 2].real,
                                         eigv_norm[i, 2].imag))
        print('     (%7.6f + %7.6fj)' % (eigv_norm[i, 3].real,
                                         eigv_norm[i, 3].imag))
        print('     (%7.6f + %7.6fj)' % (eigv_norm[i, 4].real,
                                         eigv_norm[i, 4].imag))
        print('     (%7.6f + %7.6fj)\n' % (eigv_norm[i, 5].real,
                                           eigv_norm[i, 5].imag))
    rollamp = np.absolute(eigv_norm[roll, :])
    spiralamp = np.absolute(eigv_norm[spiral, :])
    dutchamp = np.absolute(eigv_norm[dutch[0], :])
    dutchphase = np.arctan2(eigv_norm[dutch[0], :].imag, eigv_norm[dutch[0], :].real)*180/np.pi
    print('Eigenvector\tRoll Mode\tSpiral Mode\t\tDutch Roll')
    print('Component\tAmplitude\tAmplitude\tAmplitude\tPhase')
    print('----------------------------------------------------------------------')
    symbols = [u"\u0392", "p", "r", u"\u03BEy", u"\u03D5", u"\u03C8"]
    for i in range(len(rollamp)):
        print('    ' + u"\u0394" + '%s\t\t%7.6f\t%7.6f\t%7.6f     %4.2f\u00B0' % (symbols[i],
                                                                                    rollamp[i],
                                                                                    spiralamp[i],
                                                                                    dutchamp[i],
                                                                                    dutchphase[i]))

print("\n8.3 (1.0 pts)")
print("-------------------------")
HW8_3()
print("\n8.4 (3.0 pts)")
print("-------------------------")
HW8_4()
print("\n9.1 (3.0 pts)")
print("-------------------------")
HW9_1()
print("\nExample 9.2.1 (1.0 pts)")
print("-------------------------")
HWEx9_2_1()


#8.3 (/1) :
#8.4 (/3) :
#8.40 (/2) :
#9.1 (/3) :
#Ex. 9.2.1 (/1) :
#Total (/10) :
