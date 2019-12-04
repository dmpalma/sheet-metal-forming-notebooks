#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Domingo Morales Palma <dmpalma@us.es>
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches

def plot_bend_element(curvature=0.001):
    rho = 1/curvature
    l = 10
    theta = l/rho/2
    thetaG = np.rad2deg(l/rho/2)
    width = rho
    height = rho
    
    fig, ax = plt.subplots()
    arc = patches.Arc([0,-height/2], width, height, angle=90, theta1=-thetaG, theta2=thetaG)
    ax.add_patch(arc)
    ax.plot([0, rho/2*math.sin(-theta)], [-height/2, -height/2+rho/2*math.cos(-theta)], 'k:')
    ax.plot([0, rho/2*math.sin(theta)], [-height/2, -height/2+rho/2*math.cos(theta)], 'k:')
    ax.text(0, -0.5, r'$\theta = %.0f^{\circ}$' % thetaG, horizontalalignment='center')
    ax.text(0, -1, r'$\rho = %.0f$ mm' % rho, horizontalalignment='center')
    ax.text(0, -1.5, r'$1/\rho = %.3f$ mm$^{-1}$' % (1/rho), horizontalalignment='center')
    ax.axis([-5, 5, -2, 1])
    ax.set_axis_off()   
    ax.set_aspect("equal")
    plt.show()
    
def constants_plane_strain(E, nu, Y):
    Ep = E/(1-nu**2)
    Yp = Y*math.sqrt(3)/2
    return Ep, Yp

def e1(y, rho):
    return y/rho

def s1(y, rho, Ep, Yp):
    signo = math.copysign(1,e1(y,rho))
    s1e = Ep*e1(y,rho)
    return s1e if abs(e1(y,rho))<Yp/Ep else Yp*signo

def bending_char(t, Ep, Yp):
    rhoe = Ep*t/(2*Yp)
    Me = Yp*t**2/6
    Mp = 1.5*Me
    return rhoe, Me, Mp

def M(rho, rhoe, Me):
    return Me*rhoe/rho if rho>rhoe else Me*(3-(rho/rhoe)**2)/2

def plot_moment_curvature(rhoe, Me, Mp):
    x = np.arange(1/rhoe, 0.01, 1e-4)
    y = [M(1/i, rhoe, Me) for i in x]
    
    fig, ax = plt.subplots()
    ax.axhline(y=Me, color='k', ls=':', lw=0.5)
    ax.axhline(y=Mp, color='k', ls=':', lw=0.5)
    ax.axvline(x=1/rhoe, color='k', ls=':', lw=0.5)
    ax.plot([0, 1/rhoe], [0, Me], 'b-')
    ax.plot(x, y, 'b-')
    ax.text(0, Me, r'$M_e$')
    ax.text(0, Mp, r'$M_p$')
    ax.text(1/rhoe, 1, r'$(1/\rho)_e$')
    ax.text(0.004, Mp/2, r'$(1/\rho)_e=%.6f$ mm$^{-1}$' % (1/rhoe))
    ax.text(0.004, Mp/2-2.5, r'$\rho_e=%.0f$ mm' % (rhoe))
    ax.text(0.004, Mp/2-5, r'$M_e=%.1f$ Nm/m' % (Me))
    ax.text(0.004, Mp/2-7.5, r'$M_p=%.1f$ Nm/m' % (Mp))
    ax.axis([0, 0.01, 0, 1.1*Mp])
    ax.set_xlabel(r'Sheet curvature, $1/\rho$')
    ax.set_ylabel(r'Bending moment, $M$')       
    plt.show()

def plot_e1(rho, t):
    x = np.arange(-t/2, t/2, 0.01)
    y = [e1(i, rho) for i in x]
    
    fig, ax = plt.subplots()
    ax.axvline(x=0, color='k', ls=':', lw=0.5)
    ax.plot(y, x, 'b-')
    ax.axis([-0.15, 0.15, -t/2, t/2])
    ax.set_xlabel(r'Major strain, $\varepsilon_1$')
    ax.set_ylabel(r'Thickness, $t$')       
    plt.show()

def plot_s1(rho, t, Ep, Yp):
    x = np.arange(-t/2, t/2, 0.001)
    y = [s1(i, rho, Ep, Yp) for i in x]
    
    fig, ax = plt.subplots()
    ax.axvline(x=0, color='k', ls=':', lw=0.5)
    ax.plot(y, x, 'b-')
    ax.axis([-Yp-10, Yp+10, -t/2, t/2])
    ax.set_xlabel(r'Major stress, $\sigma_1$')
    ax.set_ylabel(r'Thickness, $t$')       
    plt.show()

def plot_bending1(rhoe, Me, Mp, curvature=0.001):
    rho = 1/curvature
    M1 = M(rho, rhoe, Me)
    l = 10
    theta = l/rho/2
    thetaG = np.rad2deg(l/rho/2)
    width = rho
    height = rho
    
    x = np.arange(1/rhoe, 0.25, 1e-4)
    y = [M(1/i, rhoe, Me) for i in x]
    
    fig, ax = plt.subplots(2, 1)
    arc = patches.Arc([0,-height/2], width, height, angle=90, theta1=-thetaG, theta2=thetaG)
    ax[0].add_patch(arc)
    ax[0].plot([0, rho/2*math.sin(-theta)], [-height/2, -height/2+rho/2*math.cos(-theta)], 'k:')
    ax[0].plot([0, rho/2*math.sin(theta)], [-height/2, -height/2+rho/2*math.cos(theta)], 'k:')
    ax[0].text(0, -0.5, r'$\theta = %.0f^{\circ}$' % thetaG, horizontalalignment='center')
    ax[0].text(0, -1, r'$\rho = %.0f$ mm' % rho, horizontalalignment='center')
    ax[0].text(0, -1.5, r'$1/\rho = %.3f$ mm$^{-1}$' % (1/rho), horizontalalignment='center')
    ax[0].axis([-5, 5, -2, 1])
    ax[0].set_axis_off()   
    ax[0].set_aspect("equal")

    ax[1].plot([0, 1/rhoe], [0, Me], 'b-')
    ax[1].plot(x, y, 'b-')
    ax[1].plot([curvature, curvature], [0, M1], 'r-')
    text = 'fully elastic' if M1 < Me else 'elastic+plastic'
    ax[1].text(0.08, Mp/2, r'$M=%.1f$ Nm/m (%s)' % (M1, text))
    ax[1].axis([0, 0.25, 0, 1.1*Mp])
    ax[1].set_xlabel(r'Sheet curvature, $1/\rho$')
    ax[1].set_ylabel(r'Bending moment, $M$')       
    plt.show()

def plot_bending(t, Ep, Yp, rhoe, Me, Mp, curvature=0.001):
    rho = 1/curvature
    M1 = M(rho, rhoe, Me)
    l = 10
    theta = l/rho/2
    thetaG = np.rad2deg(l/rho/2)
    width = rho
    height = rho
    
    ax1 = plt.subplot2grid((2, 3), (0, 0))
    ax2 = plt.subplot2grid((2, 3), (0, 1))
    ax3 = plt.subplot2grid((2, 3), (0, 2))
    ax4 = plt.subplot2grid((2, 3), (1, 0), colspan=3)

    arc = patches.Arc([0,-height/2], width, height, angle=90, theta1=-thetaG, theta2=thetaG)
    ax1.add_patch(arc)
    ax1.plot([0, rho/2*math.sin(-theta)], [-height/2, -height/2+rho/2*math.cos(-theta)], 'k:')
    ax1.plot([0, rho/2*math.sin(theta)], [-height/2, -height/2+rho/2*math.cos(theta)], 'k:')
#==============================================================================
#     ax1.text(0, -0.5, r'$\theta = %.0f^{\circ}$' % thetaG, horizontalalignment='center')
#     ax1.text(0, -1, r'$\rho = %.0f$ mm' % rho, horizontalalignment='center')
#     ax1.text(0, -1.5, r'$1/\rho = %.3f$ mm$^{-1}$' % (1/rho), horizontalalignment='center')
#==============================================================================
    ax1.axis([-5, 5, -2, 1])
    ax1.set_axis_off()   
    ax1.set_aspect("equal")

    x = np.arange(-t/2, t/2, 0.01)
    y = [e1(i, rho) for i in x]
    ax2.axvline(x=0, color='k', ls=':', lw=0.5)
    ax2.plot(y, x, 'r-')
    ax2.axis([-0.15, 0.15, -t/2, t/2])
    ax2.set_xlabel(r'Major strain, $\varepsilon_1$')
    ax2.set_ylabel(r'Thickness, $t$')       

    x = np.arange(-t/2, t/2, 0.001)
    y = [s1(i, rho, Ep, Yp) for i in x]
    ax3.axvline(x=0, color='k', ls=':', lw=0.5)
    ax3.plot(y, x, 'r-')
    ax3.axis([-Yp-10, Yp+10, -t/2, t/2])
    ax3.set_xlabel(r'Major stress, $\sigma_1$')
    ax3.set_ylabel(r'Thickness, $t$')
    
    x = np.arange(1/rhoe, 0.25, 1e-4)
    y = [M(1/i, rhoe, Me) for i in x]
    ax4.plot([0, 1/rhoe], [0, Me], 'b-')
    ax4.plot(x, y, 'b-')
    ax4.plot([curvature, curvature], [0, M1], 'r-')
    text = 'fully elastic' if M1 < Me else 'elastic+plastic'
    ax4.text(0.08, Mp/2, r'$M=%.1f$ Nm/m (%s)' % (M1, text))
    ax4.axis([0, 0.25, 0, 1.1*Mp])
    ax4.set_xlabel(r'Sheet curvature, $1/\rho$')
    ax4.set_ylabel(r'Bending moment, $M$')
    
    #plt.show()
    plt.tight_layout()
    
if __name__ == "__main__":
    t = 1.2
    E = 210e3
    nu = 0.3
    Y = 100
    
    Ep, Yp = constants_plane_strain(E, nu, Y)
    print('Material constants in plane strain: Ep = %.1f GPa, Yp = %.1f MPa' %(Ep/1e3, Yp))
    
    rhoe, Me, Mp = bending_char(t, Ep, Yp)
    print('Limiting elastic curvature: (1/rho)e = %.6f mm-1 --> radius = %.0f mm' % (1/rhoe, rhoe))
    print('Limiting elastic moment: Me = %.1f Nm/m' % (Me))
    print('Fully plastic moment: Mp = %.1f Nm/m' % (Mp))
    
#    plot_moment_curvature(rhoe, Me, Mp)    
#    plot_e1(4.5, t)
#    plot_s1(104.5, t, Ep, Yp)
#    plot_bending1(rhoe, Me, Mp, 0.0495)
    plot_bending(Ep, Yp, rhoe, Me, Mp, 0.0495)
