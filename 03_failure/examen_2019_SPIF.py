#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AuAuthor: Domingo Morales Palma <dmpalma@us.es>
"""

import math
import matplotlib.pyplot as plt

def plot_stress(σeff, σ1, σ2, ε1, ε2):
    alpha = σ2/σ1
    beta = ε2/ε1
    
    s1 = lambda alpha: σeff/math.sqrt(1-alpha+alpha**2)
    a0 = [200*i/1000-100 for i in range(1000)]
    yM = [s1(i) for i in a0] + [-s1(i) for i in a0]
    xM = [j*i for i,j in zip(yM,a0)] + [-j*i for i,j in zip(yM,a0)]
    
    fig, ax = plt.subplots(figsize=(6,6))
    ax.axvline(x=0, color='k', linewidth=0.2)
    ax.axhline(y=0, color='k', linewidth=0.2)
    ax.plot(xM,yM, 'k-', label='Yield surface (Mises)')
    ax.plot([0,σ2], [0,σ1], 'k:', label=r'$\alpha=%.3f$' % alpha)
    ax.plot([σ2, σ2+150*ε2], [σ1, σ1+150*ε1], color='b', label=r'$d\varepsilon_{ij}^p \quad (\beta=%.3f)$' % beta)
    ax.axis([-2*σeff, 2*σeff, -2*σeff, 2*σeff])
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\sigma_2$')
    ax.set_ylabel(r'$\sigma_1$')
    plt.legend()
    plt.show()
    
def plot_strain(ε1, ε2):
    beta = ε2/ε1
    ε3 = -(ε1+ε2)
    
    e1c = lambda b: -ε3/(1+b)
    e2c = lambda b: b*e1c(b)
    
    e2a, e1a = e2c(-0.5), e1c(-0.5)
    e2b, e1b = e2c(1), e1c(1)
    
    fig, ax = plt.subplots(figsize=(6,6))
    ax.axvline(x=0, color='k', lw=0.2)
    ax.plot([0,-1], [0,2], 'k', lw=0.2)
    ax.plot([0,1], [0,1], 'k', lw=0.2)
    ax.plot([e2a,e2b], [e1a,e1b], 'r-', label='Fracture limit')
    ax.plot([0,ε2], [0,ε1], 'b:', label=r'$\beta=%.3f$' % beta)
    ax.axis([-1, 1, 0, 2])
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\varepsilon_2$')
    ax.set_ylabel(r'$\varepsilon_1$')
    plt.legend()
    plt.show()
    
if __name__ == "__main__":
    plot_stress(308.1, 354.7, 200.9, 1.030, 0.095)
    plot_strain(1.030, 0.095)
