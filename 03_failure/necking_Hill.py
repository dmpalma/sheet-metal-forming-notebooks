#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AuAuthor: Domingo Morales Palma <dmpalma@us.es>
"""

from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

def e1c(beta, n, e0=0):
    return n/(1+beta) - e0*sqrt(3)/2*sqrt(1+beta+beta**2)


def plot_Hill(n, e0list=[]):
    b = np.linspace(-0.99, 1, 100)
    e10 = [e1c(beta, n) for beta in b]
    e20 = [beta*e1c(beta, n) for beta in b]
    
    fig, ax = plt.subplots(figsize=(6,6))
    ax.axvline(x=0, color='k', lw=0.2)
    ax.plot((0,1), (0,1), 'k', lw=0.2)
    ax.plot((0,-0.5), (0,1), 'k', lw=0.2)
    ax.plot(e20, e10, label=r'$\varepsilon_0 = 0$')
    
    for e0 in e0list:
        e11 = [e1c(beta, n, e0) for beta in b]
        e21 = [beta*e1c(beta, n, e0) for beta in b]
        
        ax.plot(e21, e11, label=r'$\varepsilon_0 = %s$' % e0)
    
    ax.axis([-0.25, 0.25, 0, 0.5])
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\varepsilon_2$')
    ax.set_ylabel(r'$\varepsilon_1$')
    if not e0list:
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
