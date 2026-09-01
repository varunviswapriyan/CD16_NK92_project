import os, sys, shutil
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import bionetgen

#from pyswarm import pso
from pyswarms.single.global_best import GlobalBestPSO
from calculate_pZAP import run_bionetgen
from interpolation import spline
from ca_ODE import calcium, solve
from interpolation_exp_data import exp_idata

from data_slice import data_after_tstart





def calculate_residue(n,tnew,PZAP,PSYK,exp_time,exp_data,C1,C2,g,N,dir_name,k3):
    
    print("C1 = "+str(C1))
    print("C2 = "+str(C2))
    print("g = "+str(g))
#3. Take the experimental data and interpolate at the theoretical points
    idata,tnew=exp_idata(exp_time, exp_data, tnew)
    idata=np.asarray(idata)
    tnew=np.asarray(tnew)
    #print(idata)
    
    
#4. Make TT, MODEL, EXPT after tstart      
    TT,EXPT,count, CA0 = data_after_tstart(tnew,idata)
    # print(count)
    # print(CA0)
    # print(len(tnew))
    # print(len(idata))

    
#5. Ca Signal from the ODE model uisng the PZAP signal   
    y0z = [CA0, 1]
    ca_PZAP,h_PZAP=solve(tnew,N,y0z,PZAP,C1,C2,g,k3)
    ca_PZAP=np.asarray(ca_PZAP)
    h_PZAP=np.asarray(h_PZAP)
    #plt.plot(tnew,ca,'m') #model generated
    print(f'ca0 in main',ca_PZAP[0])
    # plt.show()
    print(ca_PZAP)

#5. Ca Signal from the ODE model uisng the PSYK signal  

    # y0s = [CA0/2, 1]
    # ca_PSYK,h_PSYK=solve(tnew,N,y0s,PSYK,C1,C2,g)
    # ca_PSYK=np.asarray(ca_PSYK)
    # h_PSYK=np.asarray(h_PSYK)
    # #plt.plot(tnew,ca,'m') #model generated
    # print(f'ca0 in main',ca_PSYK[0])
    # # plt.show()
    # print(ca_PSYK)
    
    #ca=ca_PZAP+ca_PSYK
    ca=ca_PZAP


#6. fit starts from tstart, so, MODEL and EXPT data containts data after tstart    
    MODEL=[]
    for i in range (count,len(tnew),1): 
        MODEL.append(ca[i])

    TT=np.asarray(TT)
    EXPT=np.asarray(EXPT)
    MODEL=np.asarray(MODEL)


#7. Finding objective function only at t>25 sec
    ca_diff=(idata-ca)
    OBJ=(EXPT-MODEL)
    
    #print(ca_diff)
    SSR=np.sum(np.square(OBJ))
    #print(SSR)

    
    return SSR,TT,EXPT,MODEL    