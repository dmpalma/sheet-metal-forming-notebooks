#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Domingo Morales Palma <dmpalma@us.es>
"""

import math
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt

def sij2array(sx, sy, sz, sxy, sxz, syz):
    a = np.array([[sx, sxy, sxz],
                 [sxy, sy, syz],
                 [sxy, syz, sz]])
    return a

def I1(sx, sy, sz, sxy, sxz, syz):
    return sx + sy + sz

def I2(sx, sy, sz, sxy, sxz, syz):
    return sxy**2 + sxz**2 + syz**2 - sy*sz - sz*sx - sx*sy

def I3(sx, sy, sz, sxy, sxz, syz):
    return sx*sy*sz +2*syz*sxz*sxy - sx*syz**2 - sy*sxz**2 - sz*sxy**2
    
def calculate_invariants(sx, sy, sz, sxy, sxz, syz):
    i1 = I1(sx, sy, sz, sxy, sxz, syz)
    i2 = I2(sx, sy, sz, sxy, sxz, syz)
    i3 = I3(sx, sy, sz, sxy, sxz, syz)
    return (i1, i2, i3)

def principal_stresses(sx, sy, sz, sxy, sxz, syz):
    # method 1: computing the eigenvalues
    #a = sij2array(sx, sy, sz, sxy, sxz, syz)
    #ps = linalg.eigvals(a)
    # method 2: using the invariants to compute the cubic equation 
    i1, i2, i3 = calculate_invariants(sx, sy, sz, sxy, sxz, syz)
    ps = np.roots((1, -i1, -i2, -i3))
    return ps
    



def print_tensor(sx, sy, sz, sxy=0, sxz=0, syz=0):
    a = sij2array(sx, sy, sz, sxy, sxz, syz)
    print(a)

def plot_Mhor_circles(s1, s2, s3, sx=0, sy=0, sxy=0):
    fig, ax = plt.subplots()
    ax.axhline(y=0, color='k')
    ax.axvline(x=0, color='k')
    circle1 = plt.Circle(((s1+s2)/2, 0), abs(s1-s2)/2, fill=False)
    circle2 = plt.Circle(((s2+s3)/2, 0), abs(s2-s3)/2, fill=False)
    circle3 = plt.Circle(((s3+s1)/2, 0), abs(s3-s1)/2, fill=False)
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.add_patch(circle3)

    plt.text(s1, 0, r'$\sigma_1$', verticalalignment='top')
    plt.text(s2, 0, r'$\sigma_2$', verticalalignment='top')
    plt.text(s3, 0, r'$\sigma_3$', verticalalignment='top')
    text = r'$\sigma_1 = %.2f$ MPa, $\sigma_2 = %.2f$ MPa, $\sigma_3=%.2f$ MPa' % (s1, s2, s3)
    plt.text(0.5, 1, text, horizontalalignment='center', verticalalignment='top', transform=ax.transAxes)
    
    if not (sx==0 and sy==0 and sxy==0):
        ax.plot((sx, sy), (sxy, -sxy), 'ro')
        ax.plot((sx, sy), (sxy, -sxy), 'r--')
        text = r'$\sigma_x = %.2f$ MPa, $\sigma_y = %.2f$ MPa, $\tau_{xy}=%.2f$ MPa' % (sx, sy, sxy)
        plt.text(0.5, 0, text, color='r', horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes)
    ax.set_aspect('equal', adjustable='datalim')
    ax.plot()
    plt.show()

def mises(sx, sy, sz, sxy=0, sxz=0, syz=0):
    return math.sqrt(1/2*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+6*(sxy**2+sxz**2+syz**2)))

#def mises(s1, s2, s3):
#    return math.sqrt(1/2*((s1-s2)**2+(s2-s3)**2+(s3-s1)**2))

def tresca(s1, s2, s3):
    return max(s1, s2, s3) - min(s1, s2, s3)

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
        ax.plot([0, py], [0, px], 'r--', label=r'Stress path')
        ax.plot(py, px, 'ro')
        ax.annotate(r'(%0.1f, %0.1f) MPa' % (py,px), [py,px+20])
#        ax.annotate(r'$\alpha=%0.3f$' % alpha, [py/2,px/2])
    
    ax.set_xlabel(r'$\sigma_2$')
    ax.set_ylabel(r'$\sigma_1$')
    ax.set_aspect('equal', adjustable='datalim')
#    plt.legend()
    plt.show()


def plot_tresca_mises(Y, px, py):
    alpha = py/px
    
    tresca = ((Y, 0), (Y, Y), (0, Y), (-Y, 0), (-Y, -Y), (0, -Y), (Y, 0))
    
    a = [200*i/1000-100 for i in range(1000)]
    y0 = [Y/(math.sqrt(1-i+i**2)) for i in a]
    x0 = [j*i for i,j in zip(y0,a)]
    y1 = [-Y/(math.sqrt(1-i+i**2)) for i in a]
    x1 = [j*i for i,j in zip(y1,a)]
    
    fig, ax = plt.subplots()
    ax.axvline(x=0, color='k', linewidth=0.2)
    ax.axhline(y=0, color='k', linewidth=0.2)
    
    ax.plot(*zip(*tresca), 'g-', label=r'Tresca')
    ax.plot(x0,y0, 'b-', label=r'Mises')
    ax.plot(x1,y1, 'b-')
    ax.plot([0, py], [0, px], color='r', label=r'Stress path ($\alpha=\sigma_2/\sigma_1$)')
    
    ax.set_xlabel(r'$\sigma_2$')
    ax.set_ylabel(r'$\sigma_1$')
    ax.annotate(r'(%0.1f, %0.1f) MPa' % (py,px), [py,px+20])
    ax.annotate(r'$\alpha=%0.3f$' % alpha, [py/2,px/2])
    ax.set_aspect('equal', adjustable='datalim')
    plt.legend(title=r'$Y=%s$ MPa' % Y)
    plt.show()

if __name__ == "__main__":
    sx, sy, sz, sxy, sxz, syz = 0, 180, 75, 134, 0, 0
    print_tensor(sx, sy, sz, sxy, sxz, syz)
    w, v = eigen(sx, sy, sz, sxy, sxz, syz)
    print('E-value:', w)
    print('E-vector:\n', v)
    s1, s2, s3 = principal_stresses(sx, sy, sz, sxy, sxz, syz)
    print(s1, s2, s3)

