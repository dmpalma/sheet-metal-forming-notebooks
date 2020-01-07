#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Domingo Morales Palma <dmpalma@us.es>
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches, gridspec
from scipy.optimize import fsolve

def plot_stretching(R, TL, CL, mu, t0, K, n, angle):
    if angle == 0:
        angle = 0.001
    ca = cos(radians(angle))
    sa = sin(radians(angle))
    ta = tan(radians(angle))
    s = R*(1-ca)-ta*(R*sa-TL/2)
    xP = 0
    yP = s-R
    xA = R*sa
    yA = s-R*(1-ca)
    xM = TL/2
    yM = 0
    xB = xM - CL*ca
    yB = yM + CL*sa
    sOA = R*radians(angle)
    sAB = ((xA-xB)**2+(yA-yB)**2)**0.5
    e1pro = log((sOA+sAB)/(TL/2-CL))
    
    def equations(p):
        e1O, e1A = p
        eq1 = e1pro - ((e1O+e1A)/2*sOA+e1A*sAB)/(sOA+sAB)
        eq2 = (e1A/e1O)**n*exp(e1O-e1A) - exp(mu*radians(angle))
        return (eq1, eq2)
    if angle == 0.001:
        e1O, e1A = 0, 0
    else:
        e1O, e1A = fsolve(equations, (e1pro/2, e1pro*2))
    tO = t0*exp(-e1O)
    tA = t0*exp(-e1A)

    Kp = 2*K/sqrt(3) # plane strain
    def funcT1(e1):
        return Kp*e1**n*t0*exp(-e1)
    
    T1O = funcT1(e1O)
    T1A = funcT1(e1A)
    p = T1O/R
    F = 2*T1A*sa

    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1.5, 1])
    ax1 = plt.subplot(gs[:, 0])
    ax2 = plt.subplot(gs[0, 1])
    ax3 = plt.subplot(gs[1, 1])

#    fig, ax = plt.subplots(figsize=(14,7))
    circle = plt.Circle((0,s-R), R, ec=None, color='#cccccc')
    ax1.add_patch(circle)
    ax1.add_patch(patches.Arc((xP,yP), 2*R, 2*R, angle = 90-angle,
             theta1=0, theta2=angle, ec='b', lw=3))
    ax1.plot([xA, xB], [yA, yB], 'b-', lw=3)
    ax1.plot([xB, xM], [yB, yM], color='#cccccc', lw=20)
    note = r'Punch stroke $s = %.1f$ mm' % s
    note += '\n'
    note += r'$\theta = %.1f \degree$' % angle
    note += '\n'
    note += r'Length OA = $%.1f$ mm' % sOA
    note += '\n'
    note += r'Length AB = $%.1f$ mm' % sAB
    note += '\n'
    note += r'Total length OAB = $%.1f$ mm' % (sOA+sAB)
    note += '\n'
    note += r'Average strain $\varepsilon_1 = %.3f$' % e1pro
    note += '\n'
    note += r'Strain in O: $\varepsilon_1 = %.3f$' % e1O
    note += '\n'
    note += r'Strain in AB: $\varepsilon_1 = %.3f$' % e1A
    note += '\n'
    note += r'Thickness in O: $t = %.3f$ mm' % tO
    note += '\n'
    note += r'Thickness in AB: $t = %.3f$ mm' % tA
    note += '\n'
    note += r'Force in O: $T_1 = %.1f$ kN/m' % T1O
    note += '\n'
    note += r'Force in AB: $T_1 = %.1f$ kN/m' % T1A
    note += '\n'
    note += r'Punch pressure: $p = %.1f$ MPa' % p
    note += '\n'
    note += r'Punch force: $F = %.1f$ kN/m' % F
    ax1.annotate(note, xy=(TL/2,TL), ha='right', va='top')
    ax1.annotate('O', xy=(0,s), color='r')
    ax1.annotate('A', xy=(xA,yA), color='r')
    ax1.annotate('B', xy=(xB,yB), color='r')
    ax1.axis([0, TL/2, 0, TL])
    ax1.set_aspect("equal")
    
    # Strain and thickness
    xmax = 1.5*(TL/2 - CL)
    x = (0, sOA, sOA+sAB)
    y = (e1O, e1A, e1A)
    ax2.axvline(x=sOA, color='k', ls=':', lw=0.5)
    ax2.axvline(x=sOA+sAB, color='k', ls=':', lw=0.5)
    ax2.axhline(y=n, color='b', ls='--', lw=0.5)
    ax2.annotate('Limit strain', xy=(xmax,n), ha='right', color='b')
    ax2.plot(x, y, 'b-')
    ax2.set_ylim(0, n+0.1)
    ax2.set_xlabel('Length along the sheet')
    ax2.set_ylabel(r'Strain, $\varepsilon_1$', color='b')
    ax2p = ax2.twinx()
    y = (t0*exp(-e1O), t0*exp(-e1A), t0*exp(-e1A))
    ax2p.plot(x, y, 'c--')
    ax2p.set_xlim(0, xmax)
    ax2p.set_ylim(0, t0)
    ax2p.set_ylabel(r'Sheet thickness, $t$ (mm)', color='c')
    label = ['O', 'A', 'B']
    [ax2p.annotate(xy=[i, 0.01], s=j) for i, j in zip(x, label)]
    
    # Tension and pressure
    x = (0, sOA, sOA, sOA+sAB, sOA+sAB)
    y = [T1O, T1A, T1A, T1A, T1A]
    ax3.axvline(x=sOA, color='k', ls=':', lw=0.5)
    ax3.axvline(x=sOA+sAB, color='k', ls=':', lw=0.5)
    ax3.plot(x, y, 'm-')
    ax3.set_ylim(0, Kp)
    ax3.set_xlabel('Length along the sheet')
    ax3.set_ylabel(r'Tension, $T_1$ (kN/m)', color='m')
    ax3p = ax3.twinx()
    y = [T1O/R, T1A/R, 0, 0, 0]
    ax3p.plot(x, y, 'r--')
    ax3p.set_xlim(0, xmax)
    ax3p.set_ylim(0, 2*Kp/R)
    ax3p.set_ylabel(r'Punch pressure, $p$ (MPa)', color='r')
    
    plt.tight_layout()
    
if __name__ == "__main__":
    plot_stretching(R=1100, TL=3000, CL=300, mu=0.1, t0=1.2, K=810, n=0.24, angle=38)
