#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Python source code
"""

import math
import numpy as np
import matplotlib.pyplot as plt


def bending_variables(E, nu, S):
    Ep = E/(1-nu**2)
    Sp = S*2/math.sqrt(3)
    e1 = lambda y, rho: y/rho
    s1 = lambda y, rho: Ep*e1(y,rho)*math.copysign(1,e1(y,rho)) if abs(e1(y,rho))<Sp/Ep else Sp*math.copysign(1,e1(y,rho))
    Me = Sp*t**2/6
    Mp = 1.5*Me
    rhoe = Ep*t/(2*Sp)
    m = lambda rho: rho/rhoe
    M = lambda rho: Me if rho>rhoe else Me*(3-m(rho)**2)/2


def plot_moment_vs_curvature(rhoe, Me, Mp):
    x = np.arange(1/rhoe, 0.01, 1e-4)
    y = [M(1/i) for i in x]
    
    fig, ax = plt.subplots()
    ax.axhline(y=Me, color='k', ls=':', lw=0.5)
    ax.axhline(y=Mp, color='k', ls=':', lw=0.5)
    #ax.axvline(x=0, color='k')
    ax.plot([0, 1/rhoe], [0, Me], 'b-')
    ax.plot(x, y, 'b-')
    ax.axis([0, 0.01, 0, 1.1*Mp])
    ax.set_xlabel(r'Sheet curvature, $1/\rho$')
    ax.set_ylabel(r'Bending moment, $M$')
       
    plt.show()
