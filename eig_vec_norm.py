#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 13:23:26 2019

@author: christian
"""

import numpy as np


def normalize_evec(evec):
    evec_T = np.copy(evec.transpose())
    evec_norm = np.zeros_like(evec)
    for i in range(len(evec)):
        mag_pre = evec_T[i, :].real**2 + evec_T[i, :].imag**2
        big_loc = np.argmax(np.sqrt(mag_pre))
        conj = np.conjugate(evec_T[i, big_loc])
        evec_T[i, :] *= conj
        mag_post = evec_T[i, :].real**2 + evec_T[i, :].imag**2
        norm = np.sum(mag_post)
        evec_norm[i, :] = evec_T[i, :]/np.sqrt(norm)
    return evec_norm
