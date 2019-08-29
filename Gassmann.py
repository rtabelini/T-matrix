# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 17:36:29 2013

@author: fsantos
"""
import numpy as np

def BiotKs(Kd, Km, Kf, phi):
    gamma = 1.0 - phi - Kd/Km
    return Kd + (gamma + phi)**2/(gamma/Km + phi/Kf)

def Ks(Kd, Km, Kf, phi):
    gamma = 1.0 - phi - Kd/Km
    return Kd + (gamma + phi)**2/(gamma/Km + phi/Kf)

def Kd(Ks, Km, Kf, phi):
    gamma = phi*(Km/Kf - 1.0)
    return (Ks*(gamma + 1.0) - Km)/(gamma - 1.0 + Ks/Km)

def MKs(Kd, Km, Kf, phi):
    w = 0.5
    beta1 = (-183.05)/(1.+np.e**((phi+0.56468)/0.10817)) + 0.99494
    beta2 = 1.-(1.-phi)**3.8
    beta = (beta1**w)*(beta2**(1.-w))
#    beta = 1. - Kd/Km
    M = ((beta-phi)/Km + (phi/Kf))**(-1)
    return Km*(1.-beta) + (beta**2)*M

def MKd(Ks, Km, Kf, phi):
    beta = 1. - Kd/Km
    M = ((beta-phi)/Km + (phi/Kf))**(-1)
    return Km*(1.-beta) + (beta**2)*M