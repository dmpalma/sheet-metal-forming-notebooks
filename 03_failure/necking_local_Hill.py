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


def plot_Hill(n=0.2, e01=0.02, e02=0.04):
    b = np.linspace(-0.99, 1, 100)
    e10 = [e1c(beta, n) for beta in b]
    e20 = [beta*e1c(beta, n) for beta in b]
    e11 = [e1c(beta, n, e01) for beta in b]
    e21 = [beta*e1c(beta, n, e01) for beta in b]
    e12 = [e1c(beta, n, e02) for beta in b]
    e22 = [beta*e1c(beta, n, e02) for beta in b]
    
    fig, ax = plt.subplots(figsize=(6,6))
    ax.axvline(x=0, color='k', lw=0.2)
    ax.plot((0,1), (0,1), 'k', lw=0.2)
    ax.plot((0,-0.5), (0,1), 'k', lw=0.2)
    ax.plot(e20, e10, label=r'$\varepsilon_0 = 0$')
    ax.plot(e21, e11, label=r'$\varepsilon_0 = %s$' % e01)
    ax.plot(e22, e12, label=r'$\varepsilon_0 = %s$' % e02)
    ax.axis([-0.25, 0.25, 0, 0.5])
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\varepsilon_2$')
    ax.set_ylabel(r'$\varepsilon_1$')
    plt.title(r'Hill localized necking model, $\sigma_y=K(\varepsilon_0+\varepsilon)^n$')
    plt.legend(title=r'Material: $n=%s$' % n)
    #plt.savefig('necking_local_Hill.png')
    plt.show()
    
if __name__ == "__main__":
    plot_Hill()
    plot_Hill(n=0.25, e01=0.1, e02=0.2)
