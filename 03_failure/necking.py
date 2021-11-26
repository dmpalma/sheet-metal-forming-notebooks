#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AuAuthor: Domingo Morales Palma <dmpalma@us.es>
"""

from math import sqrt
import numpy as np
import matplotlib.pyplot as plt


def a(b):
    return (2*b+1)/(2+b)
    
# Swift with s=K*e^n

def e1s(b, n):
    return n*sqrt(3)/(2*sqrt(1+b+b**2)) * 4*(1-a(b)+a(b)**2)**(3/2) / ((2-a(b))**2 + (2*a(b)-1)**2*a(b))

def e2s(b, n):
    return b*e1s(b, n)

# Hill with s=K(e0+e)^n

def e1h(beta, n, e0=0):
    return n/(1+beta) - e0*sqrt(3)/2*sqrt(1+beta+beta**2)

def e2h(b, n, e0=0):
    return b*e1h(b, n, e0)



def plot_Swift(n):
    beta = np.linspace(-0.99, 1, 100)
    e10 = [e1s(b, n) for b in beta]
    e20 = [e2s(b, n) for b in beta]
    
    fig, ax = plt.subplots(figsize=(6,6))
    ax.axvline(x=0, color='k', lw=0.2)
    ax.plot((0,1), (0,1), 'k', lw=0.2)
    ax.plot((0,-0.5), (0,1), 'k', lw=0.2)
    ax.plot(e20, e10, label=r'$\varepsilon_0 = 0$')
    
    ax.axis([-0.25, 0.25, 0, 0.5])
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\varepsilon_2$')
    ax.set_ylabel(r'$\varepsilon_1$')
    text = r'$\sigma_y=K(\overline{\varepsilon})^n$'
    plt.title(r'Swift diffuse necking model, %s' % text)
    plt.text(0.5, 1, r'Material: $n=%s$' % n, horizontalalignment='center', verticalalignment='top', transform=ax.transAxes)
    #plt.legend(title=r'Material: $n=%s$' % n)
    #plt.savefig('necking_local_Swift.png')
    plt.show()

def plot_Swift_Hill(n):
    beta = np.linspace(-0.99, 1, 100)
    beta1 = np.linspace(-0.99, 0, 100)
    e10 = [e1s(b, n) for b in beta]
    e20 = [e2s(b, n) for b in beta]
    e11 = [e1h(b, n) for b in beta1]
    e21 = [e2h(b, n) for b in beta1]
    
    fig, ax = plt.subplots(figsize=(6,6))
    ax.axvline(x=0, color='k', lw=0.2)
    ax.plot((0,1), (0,1), 'k', lw=0.2)
    ax.plot((0,-0.5), (0,1), 'k', lw=0.2)
    ax.plot((0,-1), (0,1), 'k', lw=0.2)
    ax.plot((0,-2), (0,1), 'k', lw=0.2)
    ax.plot(e20, e10, label=r'Swift')
    ax.plot(e21, e11, label=r'Hill')
    
    ax.axis([-0.25, 0.25, 0, 0.5])
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\varepsilon_2$')
    ax.set_ylabel(r'$\varepsilon_1$')
    text = r'$\sigma_y=K(\overline{\varepsilon})^n$'
    plt.title(r'Necking models, with %s' % text)
    plt.legend(title=r'Material: $n=%s$' % n)
    #plt.savefig('necking_local_Swift.png')
    plt.show()

def plot_Hill(n, e0=[]):
    b = np.linspace(-0.99, 0, 100)
    e10 = [e1h(beta, n) for beta in b]
    e20 = [e2h(beta, n) for beta in b]
    
    fig, ax = plt.subplots(figsize=(6,6))
    ax.axvline(x=0, color='k', lw=0.2)
    ax.plot((0,1), (0,1), 'k', lw=0.2)
    ax.plot((0,-0.5), (0,1), 'k', lw=0.2)
    ax.plot(e20, e10, label=r'$\varepsilon_0 = 0$')
    
    for e0_ in e0:
        e11 = [e1h(beta, n, e0_) for beta in b]
        e21 = [e2h(beta, n, e0_) for beta in b]
        
        ax.plot(e21, e11, label=r'$\varepsilon_0 = %s$' % e0_)
    
    ax.axis([-0.25, 0.25, 0, 0.5])
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\varepsilon_2$')
    ax.set_ylabel(r'$\varepsilon_1$')
    if not e0:
        text = r'$\sigma_y=K(\varepsilon)^n$'
    else:
        text = r'$\sigma_y=K(\varepsilon_0+\varepsilon)^n$'
        plt.legend()
    plt.title(r'Hill localized necking model, %s' % text)
    plt.text(0.5, 1, r'Material: $n=%s$' % n, horizontalalignment='center', verticalalignment='top', transform=ax.transAxes)
    #plt.legend(title=r'Material: $n=%s$' % n)
    #plt.savefig('necking_local_Hill.png')
    plt.show()
    
if __name__ == "__main__":
    plot_Hill(0.25)
    plot_Hill(n=0.25, e0list=(0.1, 0.2, 0.3, 0.4))
