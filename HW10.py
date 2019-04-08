#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 08:31:20 2019

@author: christian
"""


import scipy.linalg
from eig_vec_norm import normalize_evec


def HW8_1():
    # Using 8.1.48
    m1 = 20.
    m2 = 20.
    c1 = 30.
    c2 = 15.
    k1 = 2.
    k2 = 100.
    A = [[-c1, 0., -(k1 + k2), k2],
         [0., -c2, k2, -k2],
         [1., 0., 0., 0.],
         [0., 1., 0., 0.]]
    B = [[m1, 0., 0., 0.],
         [0., m2, 0., 0.],
         [0., 0., 1., 0.],
         [0., 0., 0., 1.]]
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
        print('vec_norm =  (%7.6f + %7.6fj)' % (eigv_norm[i, 0].real,
                                                eigv_norm[i, 0].imag))
        print('            (%7.6f + %7.6fj)' % (eigv_norm[i, 1].real,
                                                eigv_norm[i, 1].imag))
        print('            (%7.6f + %7.6fj)' % (eigv_norm[i, 2].real,
                                                eigv_norm[i, 2].imag))
        print('            (%7.6f + %7.6fj)\n' % (eigv_norm[i, 3].real,
                                                  eigv_norm[i, 3].imag))


def HW8_2():
    # Using 8.1.48
    m1 = 20.
    m2 = 20.
    c1 = 30.
    c2 = 15.
    k1 = 0.
    k2 = 100.
    A = [[-c1, 0., -(k1 + k2), k2],
         [0., -c2, k2, -k2],
         [1., 0., 0., 0.],
         [0., 1., 0., 0.]]
    B = [[m1, 0., 0., 0.],
         [0., m2, 0., 0.],
         [0., 0., 1., 0.],
         [0., 0., 0., 1.]]
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
        print('vec_norm =  (%7.6f + %7.6fj)' % (eigv_norm[i, 0].real,
                                                eigv_norm[i, 0].imag))
        print('            (%7.6f + %7.6fj)' % (eigv_norm[i, 1].real,
                                                eigv_norm[i, 1].imag))
        print('            (%7.6f + %7.6fj)' % (eigv_norm[i, 2].real,
                                                eigv_norm[i, 2].imag))
        print('            (%7.6f + %7.6fj)\n' % (eigv_norm[i, 3].real,
                                                  eigv_norm[i, 3].imag))


print("\n8.1 (5.0 pts)")
print("-------------------------")
HW8_1()
print("\n8.2 (5.0 pts)")
print("-------------------------")
HW8_2()
