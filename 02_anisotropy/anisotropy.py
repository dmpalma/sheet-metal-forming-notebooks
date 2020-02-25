#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Domingo Morales Palma <dmpalma@us.es>
"""

import math
import matplotlib.pyplot as plt

def eff_stress_Mises(s1, alpha):
    ''' Mises effective stress in plane stress'''
    return s1*math.sqrt(1-alpha+alpha**2)

def eff_stress_Hosford(s1, alpha, r0, r90, a):
    ''' Hosford effective stress in plane stress'''
    return s1*( (r90+r0*alpha**a+r0*r90*(1-alpha)**a)/(r90*(1+r0)) )**(1/a)

def planar_anisotropy(r0, r90, r45=1):
    return (r0+r90-2*r45)/4

def normal_anisotropy(r0, r90, r45=1):
    return (r0+r90+2*r45)/4

def plot_ys(sy=300, r0=1.2, r90=1.8, a=8):
    s1 = lambda alpha: sy/eff_stress_Hosford(1, alpha, r0, r90, a)
    a0 = [200*i/1000-100 for i in range(1000)]
    yH = [s1(i) for i in a0] + [-s1(i) for i in a0]
    xH = [j*i for i,j in zip(yH,a0)] + [-j*i for i,j in zip(yH,a0)]
    
    # Hill
    s1 = lambda alpha: sy/eff_stress_Hosford(1, alpha, r0, r90, 2)
    yI = [s1(i) for i in a0] + [-s1(i) for i in a0]
    xI = [j*i for i,j in zip(yI,a0)] + [-j*i for i,j in zip(yI,a0)]
    # Mises
    s1 = lambda alpha: sy/eff_stress_Mises(1, alpha)
    yM = [s1(i) for i in a0] + [-s1(i) for i in a0]
    xM = [j*i for i,j in zip(yM,a0)] + [-j*i for i,j in zip(yM,a0)]
    
    plt.rcParams["figure.figsize"] = (6,6)
    fig, ax = plt.subplots()
    ax.axvline(x=0, color='k', linewidth=0.2)
    ax.axhline(y=0, color='k', linewidth=0.2)
    ax.plot(xM,yM, 'g-', label=r'Mises')
    ax.plot(xI,yI, 'r--', label=r'Hill')
    ax.plot(xH,yH, 'b:', label=r'Hosford')
    ax.annotate(r'$\Delta r = %0.3f$' % planar_anisotropy(r0, r90), xy=(500,700))
    ax.annotate(r'$\overline{r} = %0.3f$' % normal_anisotropy(r0, r90), xy=(500,600))
    ax.axis([-800, 800, -800, 800])
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\sigma_2$')
    ax.set_ylabel(r'$\sigma_1$')
    plt.legend(title='Yield surface', loc='lower right')
    plt.show()

if __name__ == "__main__":
    plot_ys()
    