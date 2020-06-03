# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 18:12:08 2020

@author: AX201 GMRS
"""

import numpy as np
from scipy.interpolate import interp1d

def model(z , t , M , C , K , sismo):
    """
    z would be the solution vector
    dzdt stores the system of ODEs ready for integration.
    """
    U = z[0]
    V = z[1]
    f_t  = inertial(M , sismo , t)
    dUdt = V
    dVdt = (1/M)*(f_t - C*V - K*U)
    dzdt = [dUdt , dVdt]
    
    return dzdt

def inertial(M , sismo , t ):
    """
    Computes the inertial force using interpolation
    based on the original time signal
    
    """
    ndats = len(sismo)
    t_max = (ndats-1)*0.02
    tt  = np.linspace(0 , t_max , ndats)
    accel = interp1d(tt , sismo)
    accel_new = accel(t)    
    f_t =-M*accel_new
    return f_t

def misses(F_0 , fac , lmda , x):
    """
    Constitutive law for the Misses truss. This is a piecewise
    continuous function according to fac.
    F_0 : force amplitude
    fac : number of wavelengths at which the constant slope phase starts
    lmda: Wavelength
    x   : Displacement
    """
    dx = lmda/100
    df = F_0*np.sin((2*np.pi/lmda)*dx)
    k = df/dx
    if x > -fac*lmda and x < fac*lmda:
        F = F_0*np.sin((2*np.pi/lmda)*x)
    else:
        if x<= -fac*lmda:
            F = k*(x + fac*lmda)
        else:
            F = k*(x - fac*lmda)
    return F

def dissip_device(z , t , M , C , K , sismo , m , c , k , F_0 , fac , lmda):
    """
    z would be the solution vector
    dzdt stores the system of ODEs ready for integration.
    """
    U = z[0]
    V = z[1]
    u = z[2]
    v = z[3]
#
    F_AP  = inertial(M , sismo , t)
    F_VM  = misses(F_0 , fac , lmda ,U-u)
#
    dUdt = V
    dVdt = (1/M)*(F_AP - C*V - K*U-F_VM)
    dudt = v
    dvdt = (1/m)*(F_VM - k*u - c*v)
#
    dzdt = [dUdt , dVdt , dudt , dvdt]
    
    return dzdt

def plot_misses(F_0 , lmda , ndats , fac):
    fac   = 1.0                                 #Piecewise continous at fac*lambda
    x  = np.linspace(-2*lmda , 2*lmda, ndats)   #Assumed displacements span
    n = len(x)
    FM = np.zeros(n)
#
    for i in range(n):
        FM[i] = misses(F_0 , fac , lmda , x[i])
    
    return FM , x