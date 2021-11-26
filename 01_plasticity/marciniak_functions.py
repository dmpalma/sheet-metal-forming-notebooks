#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Domingo Morales Palma <dmpalma@us.es>
"""

import math
import matplotlib.pyplot as plt


def plot_strains(e1, e2):
    beta = e2/e1
    eeff = e1*2/math.sqrt(3)*math.sqrt(1+beta+beta**2)

    # yield surface
    b = [3*i/1000-2 for i in range(1000)]
    yi = [eeff/(2/math.sqrt(3)*math.sqrt(1+i+i**2)) for i in b]
    xi = [j*i for i,j in zip(yi,b)]

    fig, ax = plt.subplots()
    ax.plot([0,-2], [0,1], 'k-', linewidth=0.2)
    ax.plot([0,0], [0,1], 'k-', linewidth=0.2)
    ax.plot([0,1], [0,1], 'k-', linewidth=0.2)
    ax.plot([0, e2], [0, e1], color='r', label='Strain path')
    ax.plot(xi,yi, 'b--', label=r'Mises yield surface at $\overline{\varepsilon}=%0.3f$' % eeff)
    ax.axis([-0.3, 0.2, 0, 0.4])
    ax.set_xlabel(r'$\varepsilon_2$')
    ax.set_ylabel(r'$\varepsilon_1$')
    ax.annotate(r'(%0.3f, %0.3f)' % (e2,e1), [e2,e1+0.01])
    ax.annotate(r'$\beta=%0.3f$' % beta, [e2/2,e1/2])
    plt.legend()
    plt.show()

def plot_stresses(s1, s2):

    alpha = s2/s1
    seff = s1*math.sqrt(1-alpha+alpha**2)

    # yield surface
    a = [200*i/1000-100 for i in range(1000)]
    yi = [seff/(math.sqrt(1-i+i**2)) for i in a]
    xi = [j*i for i,j in zip(yi,a)]

    fig, ax = plt.subplots()
    ax.axvline(x=0, color='k', linewidth=0.2)
    ax.axhline(y=0, color='k', linewidth=0.2)
    ax.plot([0, s2], [0, s1], color='r', label='Stress path')
    ax.plot(xi,yi, 'b--', label='Yield surface at $\overline{\sigma}=%0.1f$ MPa' % seff)
    ax.axis([-500, 600, -500, 600])
    ax.set_xlabel(r'$\sigma_2$')
    ax.set_ylabel(r'$\sigma_1$')
    ax.annotate(r'(%0.1f, %0.1f) MPa' % (s2,s1), [s2,s1+20])
    ax.annotate(r'$\alpha=%0.3f$' % alpha, [s2/2,s1/2])
    plt.legend(loc='lower right')
    plt.show()
    
def plot_Hollomon(K, n, eeff, seff):
    x = [i/500 for i in range(500)]
    y = [K*i**n for i in x]

    x1 = [eeff*i/500 for i in range(500)]
    y1 = [K*i**n for i in x1]

    fig, ax = plt.subplots()
    ax.plot(x, y, 'b-')
    ax.plot(x1, y1, 'r-')
    ax.axis([0, 0.4, 0, 500])
    ax.set_xlabel(r'True strain, $\varepsilon_{eff}$')
    ax.set_ylabel('True stress, $\sigma_{eff}$')
    ax.plot([eeff,eeff], [0,seff], 'r:')
    ax.annotate(r'$\varepsilon_{eff}=%0.2f$' % eeff, [eeff,10])
    ax.plot([0,eeff], [seff,seff], 'r:')
    ax.annotate(r'$\sigma_{eff}=%0.1f$ MPa' % seff, [0.01,seff])
    plt.show()
    
def plot_stress_paths(a1=-1, a2=0, a3=1/2, Y = 100):
    s1f = 300

    # yield surface
    a = [200*i/1000-100 for i in range(1000)]
    y0 = [Y/(math.sqrt(1-i+i**2)) for i in a]
    x0 = [j*i for i,j in zip(y0,a)]

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.axvline(x=0, color='k', linewidth=0.2)
    ax.axhline(y=0, color='k', linewidth=0.2)
    ax.plot(x0,y0, 'b:', label=r'Yield surface at $\overline{\sigma}=%s$ MPa' % Y)
    ax.plot([0, a1*s1f], [0, s1f], color='r', label=r'Constant thickness ($\alpha$=%s)' % a1)
    ax.plot([0, a2*s1f], [0, s1f], color='b', label=r'Uniaxial tension ($\alpha$=%s)' % a2)
    ax.plot([0, a3*s1f], [0, s1f], color='g', label=r'Plane strain ($\alpha$=%s)' % a3)
    ax.axis([-200, 200, -200, 300])
    ax.set_xlabel(r'$\sigma_2$')
    ax.set_ylabel(r'$\sigma_1$')
    plt.legend(loc='lower right')
    ax.set_aspect('equal', adjustable='datalim')
    plt.show()
    
def plot_strain_paths(b1=-1, b2=-0.5, b3=0):
    e1f = 0.2

    fig, ax = plt.subplots()
    ax.axvline(x=0, color='k', linewidth=0.2)
    ax.axhline(y=0, color='k', linewidth=0.2)
    ax.plot([0, b1*e1f], [0, e1f], color='r', label=r'Constant thickness ($\beta$=%0.1f)' % b1)
    ax.plot([0, b2*e1f], [0, e1f], color='b', label=r'Uniaxial tension ($\beta$=%0.1f)' % b2)
    ax.plot([0, b3*e1f], [0, e1f], color='g', label=r'Plane strain ($\beta$=%0.1f)' % b3)
    ax.axis([-0.3, 0.1, -0.1, 0.2])
    ax.set_xlabel(r'$\varepsilon_2$')
    ax.set_ylabel(r'$\varepsilon_1$')
    plt.legend(loc='lower left')
    ax.set_aspect('equal', adjustable='datalim')
    plt.show()

