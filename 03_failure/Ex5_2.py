#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Domingo Morales Palma <dmpalma@us.es>
"""

import math
import matplotlib.pyplot as plt


def plot_force(K, n, t0, w0A, w0B, F1, e1A):
    F1A = lambda x: K*t0*w0A * x**n * math.exp(-x)
    F1B = lambda x: K*t0*w0B * x**n * math.exp(-x)

    x = [i/1000 for i in range(1000)]
    yA = [F1A(i) for i in x]
    yB = [F1B(i) for i in x]
    
    fig, ax = plt.subplots()
    ax.plot(x, yA, 'k-', label="A")
    ax.plot(x, yB, 'k--', label="B")
    ax.plot([n, n], [0, F1], 'r:', label="Necking in B")
    ax.plot([0, n], [F1A(e1A), F1A(e1A)], 'r--', label="Maximum force in B")
    ax.plot([e1A, e1A], [0, F1A(e1A)], 'b:', label="Major strain in A")
    ax.axis([0, 0.4, F1*7/8, F1+500])
    ax.set_xlabel(r'$\varepsilon_1$')
    ax.set_ylabel(r'$F_1$ (N)')
    plt.legend()
    plt.show()
