#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 22:17:27 2022

@author: ixn004
"""

import numpy as np
from scipy.integrate import odeint

def calcium(y, t, PZAP,C1,C2,g,k3,k4):
    
    #z=602 #from uM unit to molecules/mu^3 unit
   
    be=0
    b=0.111 # maximum Ca flow from Cyto to ER  
    k1=0.7 # In Atri 1993, has unit of uM 
    k2=0.7 # In Atri 1993, has unit of uM 
    s=(k2*k2) #
    
    #K3=np.max(PZAP)
    ca, h = y
    #pZAP is in the uM unit
    
    # dcadt=(C1*h*PZAP*(((b*k1)+ca)/(k1+ca)))-(g*ca)+(be)
    # dhdt=(C2*PZAP)*((s/(s+ca ** 2))-h)
    nn=9
    #Blocked the above part on Aug 11, and changed pZAP to 
    #dcadt=(C1*h*(PZAP**nn/(PZAP**nn + k3**nn)))*(((b*k1)+ca)/(k1+ca))-(g*ca)+(be)
    #dhdt=(C2*((PZAP**nn/(PZAP**nn + k4**nn))))*((s/(s+ (ca ** 2)))-h)
    #dhdt=C2*((s/(s+ (ca ** 2)))-h)
    #k4=0.00

    #dcadt=(C1*h*((PZAP**nn/(PZAP**nn + k3**nn))+k4))*(((b*k1)+ca)/(k1+ca))-(g*ca)+(be)
    #dhdt=(C2*((PZAP**nn/(PZAP**nn + k3**nn))+k4))*((s/(s+ (ca ** 2)))-h)

    dcadt=(C1*h*((PZAP**nn/(PZAP**nn + k3**nn))+(k4*PZAP)))*(((b*k1)+ca)/(k1+ca))-(g*ca)+(be)
    dhdt=(C2*((PZAP**nn/(PZAP**nn + k3**nn))+(k4*PZAP)))*((s/(s+ (ca ** 2)))-h)

    dydt = [dcadt,dhdt]
    return dydt

def solve(tnew,N,y0,PZAP,C1,C2,g,k3,k4):   
# store solution
    #empty_like_sets_array_of_same_shape_and_size
    ca = np.empty_like(tnew)
    h = np.empty_like(tnew)
    # record initial conditions
    ca[0] = y0[0]
    h[0] = y0[1]
    #PZAP=np.ones(len(PZAP))

    # solve ODE
    for i in range(1,N):
        # span for next time step
        tspan = [tnew[i-1],tnew[i]]
        # solve for next step
        sol = odeint(calcium,y0,tspan,args=(PZAP[i],C1,C2,g,k3,k4))
        # store solution for plotting
        ca[i] = sol[1][0]
        h[i] = sol[1][1]
        # next initial condition
        y0 = sol[1]
        
    return ca,h    