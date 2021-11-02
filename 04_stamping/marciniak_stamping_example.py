#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Domingo Morales Palma <dmpalma@us.es>
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def get_vars(params, length, theta, tension):
    Rf, Rp, Rd, mu, K, n, t0 = params
    sOA, sAB, sBC, sCD, sDE, sEF = length
    T1O, T1A, T1B, T1C, T1D, T1E, T1F = tension
    thetaOA, thetaAB, thetaOB, thetaDC = theta
    
    sO = 0
    sA = sO + sOA
    sB = sA + sAB
    sC = sB + sBC
    sD = sC + sCD
    sE = sD + sDE
    sF = sE + sEF
    position = [sO, sA, sB, sC, sD, sE, sF]

    pthetaOA = [0 + thetaOA*i/9 for i in range(10)]
    pthetaAB = [thetaOA + thetaAB*i/19 for i in range(20)]
    pthetaCD = [thetaOB - thetaDC*i/19 for i in range(20)]
    psOA = [i*Rf for i in pthetaOA]
    psAB = [i*Rp+sA for i in pthetaAB]
    psCD = [i*Rd+sC for i in pthetaAB]
    psEF = [sE, sF]
    pT1OA = [T1O*math.exp(mu*i) for i in pthetaOA]
    pT1AB = [T1O*math.exp(mu*i) for i in pthetaAB]
    pT1CD = [T1D*math.exp(mu*i) for i in pthetaCD]
    pT1EF = [T1E, T1F]
    ps = psOA + psAB + psCD + psEF
    pT1 = pT1OA + pT1AB + pT1CD + pT1EF

    B = T1E/(2*mu)

    ppOA = [i/Rf for i in pT1OA]
    ppAB = [i/Rp for i in pT1AB]
    ppBC = [0, 0]
    ppCD = [i/Rd for i in pT1CD]
    ppDF = [B/sEF, B/sEF, 0]
    pp = ppOA + ppAB + ppBC + ppCD + ppDF
    ps2 = psOA + psAB + [sB, sC] + psCD + [sD, sF, sF]
    
    return position, ps, pT1, pp, ps2

def plot_T1(params, length, theta, tension):
    position, ps, pT1, pp, ps2 = get_vars(params, length, theta, tension)
    sO, sA, sB, sC, sD, sE, sF = position
    T1O, T1A, T1B, T1C, T1D, T1E, T1F = tension
    
    fig, ax = plt.subplots()
    [ax.axvline(x=i, color='grey', linestyle=':') for i in (sO, sA, sB, sC, sD)]
    ax.plot(ps, pT1, 'b-')
    ax.plot(position, tension, 'bo')
    ax.set_ylim(0, max(tension)+10)
    ax.set_xlabel('Position along the sheet (mm)')
    ax.set_ylabel(r'Tension, $T_1$', color='b')

    label = ['O', 'A', 'B', 'C', 'D', 'E', 'F']
    [plt.annotate(xy=[i, 0.01], text=j) for i, j in zip(position, label)]
    plt.show()

def plot_strain(params, length, theta, tension, strain, thickness):
    Rf, Rp, Rd, mu, K, n, t0 = params
    position, ps, pT1, pp, ps2 = get_vars(params, length, theta, tension)
    sO, sA, sB, sC, sD, sE, sF = position
    e1O, e1A, e1B, e1C, e1D, e1E, e1F = strain

    pstrain = [e1O]
    s1 = lambda e1: K*e1**n
    funcT1 = lambda x : s1(x)*t0*math.exp(-x)
    for i in range(len(pT1)-2):
        t1 = pT1[i+1]
        e1 = pstrain[i]
        funct1 = lambda x : t1 - funcT1(x)
        pstrain.append(fsolve(funct1, e1))
    pstrain.append(e1F)
    #pstrain = [fsolve(funcA, e1O) for i in pT1]

    fig, ax = plt.subplots()
    [ax.axvline(x=i, color='grey', linestyle=':') for i in (sO, sA, sB, sC, sD)]
    ax.plot(ps, pstrain, 'b--')
    ax.plot(position, strain, 'bo')
    ax.set_ylim(0, max(strain)+0.02)
    ax.set_xlabel('Position along the sheet (mm)')
    ax.set_ylabel(r'Strain, $\varepsilon_1$', color='b')

    ax2 = ax.twinx()
    ax2.plot(position, thickness, 'ro-')
    ax2.set_ylim(0, t0)
    ax2.set_ylabel(r'Thickness, $t$ (mm)', color='r')

    label = ['O', 'A', 'B', 'C', 'D', 'E', 'F']
    [plt.annotate(xy=[i, 0.01], text=j) for i, j in zip(position, label)]
    plt.show()

def plot_pressure(params, length, theta, tension, pressure):
    Rf, Rp, Rd, mu, K, n, t0 = params
    position, ps, pT1, pp, ps2 = get_vars(params, length, theta, tension)
    sO, sA, sB, sC, sD, sE, sF = position
    position2 = [sO, sA, sA, sB, sB, sC, sC, sD, sD, sE, sE, sF, sF]

    fig, ax = plt.subplots()
    [ax.axvline(x=i, color='grey', linestyle=':') for i in (sO, sA, sB, sC, sD)]
    ax.plot(ps2, pp, 'b-')
    ax.plot(position2, pressure, 'bo')
    ax.set_ylim(0, max(pressure)+2)
    ax.set_xlabel('Position along the sheet (mm)')
    ax.set_ylabel(r'Pressure, $p$', color='b')
    
    label = ['O', 'A', 'B', 'C', 'D', 'E', 'F']
    [plt.annotate(xy=[i, 0.01], text=j) for i, j in zip(position, label)]
    plt.show()


if __name__ == "__main__":
    t0 = 0.8
    thetaOB = math.pi/2
    K = 750
    n = 0.23
    s1 = lambda e1: K*e1**n   # material behaviour as a function of strain e1
    a = 330
    Rf = 2800
    Rp = 10
    Rd = 10
    sBC = 28
    sDE = 0
    sEF = 80
    e1O = 0.03
    mu = 0.1
    
    s1O = s1(e1O)
    e3O = -e1O
    tO = t0*math.exp(e3O)
    T1O = s1O*tO

    thetaOA = math.asin((a-Rp)/Rf)
    sOA = Rf*thetaOA
    thetaAB = thetaOB - thetaOA
    sAB = Rp*thetaAB
    thetaDC = thetaOB
    sCD = Rd*thetaDC
    
    T1A = T1O*math.exp(mu*thetaOA)
    T1B = T1O*math.exp(mu*3.1416/2)
    T1C = T1B
    T1D = T1C*math.exp(-mu*thetaDC)
    T1E = T1D
    T1F = 0
    
    params = [Rf, Rp, Rd, mu, K, n, t0]
    length = [sOA, sAB, sBC, sCD, sDE, sEF]
    theta = [thetaOA, thetaAB, thetaOB, thetaDC]
    tension = [T1O, T1A, T1B, T1C, T1D, T1E, T1F]

    position, ps, pT1, pp, ps2 = get_vars(params, length, theta, tension)
    sO, sA, sB, sC, sD, sE, sF = position
    plot_T1(params, length, theta, tension)
    
    B = T1E/(2*mu)
    F = 2*T1B*math.sin(thetaOB)
    
    funcT1 = lambda x : s1(x)*t0*math.exp(-x)
    funcA = lambda x : T1A - funcT1(x)
    funcB = lambda x : T1B - funcT1(x)
    funcD = lambda x : T1D - funcT1(x)
    e1A = fsolve(funcA, e1O)
    e1B = fsolve(funcB, e1A)
    e1C = e1B
    e1D = fsolve(funcD, e1C)
    e1E = e1D
    e1F = 0

    strain = [e1O, e1A, e1B, e1C, e1D, e1E, e1F]
    thickness = [t0*math.exp(-e) for e in strain]
    plot_strain(params, length, theta, tension, strain, thickness)

    pO = T1O/Rf
    pA1 = T1A/Rf
    pA2 = T1A/Rp
    pB1 = T1B/Rp
    pB2 = 0
    pC1 = 0
    pC2 = T1C/Rp
    pD1 = T1D/Rd
    pD2 = 0
    pE1 = 0
    pE2 = B/sEF
    pF1 = pE2
    pF2 = 0
    
    pressure = [pO, pA1, pA2, pB1, pB2, pC1, pC2, pD1, pD2, pE1, pE2, pF1, pF2]
    plot_pressure(params, length, theta, tension, pressure)

