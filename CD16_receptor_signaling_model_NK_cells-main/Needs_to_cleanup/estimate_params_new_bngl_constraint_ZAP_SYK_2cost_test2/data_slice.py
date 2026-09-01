#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 12:17:36 2022

@author: ixn004
"""


import numpy as np

# def data_after_tstart(tnew,idata):
#     count=0
#     for i in range(len(tnew)):
#         if (tnew[i] < 25.0):
#             count=count+1
         
#     #print(count) 
    
#     #initial point for Ca Model at t=0
#     CA0=idata[count+1] 
#     #print(CA0)
    
    
#     TT=[]
#     EXPT=[]
    
#     for i in range (count,len(tnew),1):
#         TT.append(tnew[i])
#         EXPT.append(idata[i])
#     return TT,EXPT, count, CA0    


def data_after_tstart(tnew,idata):
    count=0
    
    #initial point for Ca Model at t=0
    CA0=idata[count]
    #print(CA0)   
    TT=[]
    EXPT=[]
    
    for i in range (count,len(tnew),1):
        TT.append(tnew[i])
        EXPT.append(idata[i])
    return TT,EXPT, count, CA0   