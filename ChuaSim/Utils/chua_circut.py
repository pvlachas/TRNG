#!/usr/bin/env python
# # -*- coding: utf-8 -*-

"""Created by: Vlachas Pantelis, CSE-lab, ETH Zurich
"""
import numpy as np

"""
Dimensionless state equations, parametrization according to:
Parlitz, U. (1993). Lyapunov exponents from Chua's circuit. Journal of Circuits, Systems, and Computers, 3(02), 507-523.
"""


def chua(u, t,
         alpha = 11.852589641434262,
         beta = 14.904764875,
         gamma = 0.2977975,
         a=-1.1399612800364425,
         b=-0.7120006130847537,
         ):
    dudt = np.zeros(np.shape(u))
    x, y, z = u  
    fx = b * x + 0.5 * (a-b)*(np.abs(x + 1) - np.abs(x - 1))
    dudt[0] = alpha * (y - x - fx)
    dudt[1] = x - y + z
    dudt[2] = -beta * y - gamma * z
    return dudt


# R   = 1001
# R0  = 20
# Ga  = -1/878.1
# Gb  = 1/2339-1/878.1
# Bp  = 1
# L   = 12e-3
# C1  = 15.06e-9
# C2  = 178.5e-9


# alpha   = C2 / C1
# beta    = R*R*C2/L
# gamma   = R*R0*C2/L

# a       = R*Ga
# b       = R*Gb

# print("alpha={:}".format(alpha))
# print("beta={:}".format(beta))
# print("gamma={:}".format(gamma))
# print("a={:}".format(a))
# print("b={:}".format(b))

# alpha=11.852589641434262
# beta=14.904764875
# gamma=0.2977975
# a=-1.1399612800364425
# b=-0.7120006130847537








from scipy.integrate import odeint

def propagateChuaDynamics(u0, t_f):
    # time discretization
    t0 = 0
    dt = 1e-3
    t = np.arange(t0, t_f+dt, dt)
    # integrate ode system
    sol = odeint(chua, u0, t)
    return sol[-1]








