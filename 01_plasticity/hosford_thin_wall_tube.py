#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Domingo Morales Palma <dmpalma@us.es>
"""

import math
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
from matplotlib import patches, gridspec
from scipy.optimize import fsolve
import functions as f

         # Units:
Y = 250  # MPa
L = 2000 # mm
t = 2    # mm
D = 80   # mm
F = 8000 # N
T = 2700 # Nm

def _set_param(Y1, t1, D1, F1, T1):
    global Y, t, D, F, T
    Y = Y1
    t = t1
    D = D1
    F = F1
    T = T1

sr = 0

def st(p):
    return p*D/(2*t)
    
def sz(p):
    return p*D/(4*t) + F/(math.pi*D*t)
    
def srt(p):
    return 2*T/(math.pi*D**2*t) *1000

def compute_p(params):
    _set_param(*params)
    func = lambda p: f.mises(sr, st(p), sz(p), srt(p), 0, 0) - Y
    py, = fsolve(func, 10)
    return py

def eigen(sx, sy, sz, sxy, sxz, syz):
    a = np.array([[sx, sxy, sxz],
                 [sxy, sy, syz],
                 [sxy, syz, sz]])
    eigenValues, eigenVectors = linalg.eig(a)
    
    # sort eigenvalues and associated eigenvectors
    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    
    return eigenValues, eigenVectors

def plot_mohr(params, p):
    _set_param(*params)

    print('Stress tensor (in MPa):')
    f.print_tensor(sr, st(p), sz(p), srt(p), 0, 0)
    
    eigenValues, eigenVectors = eigen(sr, st(p), sz(p), srt(p), 0, 0)
    s1, s2, s3 = eigenValues
    print('Stress tensor in principal directions (in MPa):')
    f.print_tensor(s1, s2, s3)
    print('Principal directions:')
    print(eigenVectors)

    print("Mohr's circles (where x: radial direction, y: circumferential direction):")
    f.plot_Mhor_circles(s1, s2, s3, sr, st(p), srt(p))

def plot_mises(Y, sz=0, px=0, py=0):
    a = [200*i/1000-100 for i in range(1000)]
    y0 = [Y/(math.sqrt(1-i+i**2)) for i in a]
    x0 = [j*i for i,j in zip(y0,a)]
    y1 = [-Y/(math.sqrt(1-i+i**2)) for i in a]
    x1 = [j*i for i,j in zip(y1,a)]
    
    fig, ax = plt.subplots()
    ax.axvline(x=0, color='k', linewidth=0.2)
    ax.axhline(y=0, color='k', linewidth=0.2)    

    if not sz==0: # translation to the plane sz
        x0 = [sz + i for i in x0]
        y0 = [sz + i for i in y0]
        x1 = [sz + i for i in x1]
        y1 = [sz + i for i in y1]
        plt.text(1, 1, r'$\sigma_z=%0.2f$' % sz, color='r', horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
    
    ax.plot(x0,y0, 'b-', label=r'Mises, $Y=%s$ MPa' % Y)
    ax.plot(x1,y1, 'b-')
    plt.text(0, 1, r'Mises, $Y=%s$ MPa' % Y, color='b', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
    
    if not (px==0 and py==0):
        alpha = py/px
#        ax.plot([0, py], [0, px], 'r--', label=r'Stress path')
        ax.plot(py, px, 'ro')
        ax.annotate(r'(%0.1f, %0.1f) MPa' % (py,px), [py,px+20])
#        ax.annotate(r'$\alpha=%0.3f$' % alpha, [py/2,px/2])
    
    ax.set_xlabel(r'$\sigma_2$')
    ax.set_ylabel(r'$\sigma_1$')
    ax.set_aspect('equal', adjustable='datalim')
#    plt.legend()
    plt.show()

def plot_mohr_mises(Y, t, D, F, T, p):
    _set_param(Y, t, D, F, T)

    func = lambda p: f.mises(sr, st(p), sz(p), srt(p), 0, 0) - Y
    py, = fsolve(func, 10)

    print('Stress tensor (in MPa):')
    f.print_tensor(sr, st(p), sz(p), srt(p), 0, 0)

    eigenValues, eigenVectors = eigen(sr, st(p), sz(p), srt(p), 0, 0)
    s1, s2, s3 = eigenValues
    print('Stress tensor in principal directions (in MPa):')
    f.print_tensor(s1, s2, s3)
    print('Principal directions:')
    print(eigenVectors)

    print('Internal pressure that causes yielding: p = %.2f MPa.' % py)

    f.plot_Mhor_circles(s1, s2, s3, sr, st(p), srt(p))
    plot_mises(Y, s2, s1, s3)
    
if __name__ == "__main__":
    Y = 250
    t = 2
    D = 80
    F = 8000
    T = 2000
    p = 0
    params = (Y, t, D, F, T)
    
    plot_mohr(params, p)
    #plot_mohr_mises(Y, t, D, F, T, p)
