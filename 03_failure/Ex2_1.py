#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Domingo Morales Palma <dmpalma@us.es>
"""

import math
import matplotlib.pyplot as plt


def plot_strains(e1, e2):
    beta = e2/e1
    ey = e1*2/math.sqrt(3)*math.sqrt(1+beta+beta**2)

    # yield surface
    b = [3*i/1000-2 for i in range(1000)]
    yi = [ey/(2/math.sqrt(3)*math.sqrt(1+i+i**2)) for i in b]
    xi = [j*i for i,j in zip(yi,b)]

    fig, ax = plt.subplots()
    ax.plot([0,-2], [0,1], 'k-', linewidth=0.2)
    ax.plot([0,0], [0,1], 'k-', linewidth=0.2)
    ax.plot([0,1], [0,1], 'k-', linewidth=0.2)
    ax.plot([0, e2], [0, e1], color='r', label='Strain path')
    ax.plot(xi,yi, 'b--', label=r'Yield surface (isotropic hardening, $\varepsilon_y=%0.3f$)' % ey)
    ax.axis([-0.3, 0.2, 0, 0.4])
    ax.set_xlabel(r'$\varepsilon_2$')
    ax.set_ylabel(r'$\varepsilon_1$')
    ax.annotate(r'(%0.3f, %0.3f)' % (e2,e1), [e2,e1+0.01])
    ax.annotate(r'$\beta=%0.3f$' % beta, [e2/2,e1/2])
    plt.legend()
    plt.show()

def plot_stresses(s1, s2, s0 = 100):

    alpha = s2/s1
    sy = s1*math.sqrt(1-alpha+alpha**2)

    # yield surface
    a = [200*i/1000-100 for i in range(1000)]
    y0 = [s0/(math.sqrt(1-i+i**2)) for i in a]
    x0 = [j*i for i,j in zip(y0,a)]
    yi = [sy/(math.sqrt(1-i+i**2)) for i in a]
    xi = [j*i for i,j in zip(yi,a)]

    fig, ax = plt.subplots()
    ax.axvline(x=0, color='k', linewidth=0.2)
    ax.axhline(y=0, color='k', linewidth=0.2)
    ax.plot(x0,y0, 'b:', label=r'Initial yield surface ($\sigma_y=%s$ MPa)' % s0)
    ax.plot([0, s2], [0, s1], color='r', label='Stress path')
    ax.plot(xi,yi, 'b--', label='Yield surface (isotropic hardening, $\sigma_y=%0.1f$ MPa)' % sy)
    ax.axis([-500, 600, -500, 600])
    ax.set_xlabel(r'$\sigma_2$')
    ax.set_ylabel(r'$\sigma_1$')
    ax.annotate(r'(%0.1f, %0.1f) MPa' % (s2,s1), [s2,s1+20])
    ax.annotate(r'$\alpha=%0.3f$' % alpha, [s2/2,s1/2])
    plt.legend(loc='lower right')
    plt.show()
