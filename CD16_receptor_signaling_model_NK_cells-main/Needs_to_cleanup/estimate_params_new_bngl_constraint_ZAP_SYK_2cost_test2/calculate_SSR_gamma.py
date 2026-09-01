import os, sys, shutil
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import bionetgen

#from pyswarm import pso

from calculate_pZAP import run_bionetgen
from interpolation import spline
from ca_ODE import calcium, solve
from interpolation_exp_data import exp_idata
from data_slice import data_after_tstart



    
def gamma_PZAP(n,a0,ZAP0,SYK0,dir_name):
    
    
            
    # p=8 # which column index which obseravable in bionetgen file for total ZAP70
    # p1=13 # which column index which obseravable in bionetgen file for total SYK

    p=6 # which column index which obseravable in bionetgen file for bound ZAP70
    p1=11 # which column index which obseravable in bionetgen file for bound SYK


    N=2000 #Ca signal time points and SSR at N timepoints

    tstart=0.0 # Fit to be start from which timepoint : interpolation of pZAP70 signal starts at 0 sec

    Vc=25.0 #pZAP molecules in the simulation box of size 25 um^3
    z=602.0 #constant factor to convert from molecules/um3 to uM
    
 
    print("A0= "+str(a0))
    print("ZAP0= "+str(ZAP0))
    print("SYK0= "+str(SYK0))

    dir_name = "gamma_runs/" + str(dir_name)
    model = bionetgen.bngmodel("JJ_gamma_0_July2.bngl")

    model.parameters.adaptor0 = a0 # setting parameter A0 to a0
    model.parameters.Z0 = ZAP0 # setting parameter zeta0 to 0.0
    model.parameters.SS0 = SYK0  # setting parameter ZAP0 to 0.0

    #print(model)
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    #print model in directiry name_gamma_HPC.bngl    5 folders
    with open(f"{dir_name}/model.bngl", "w") as f:
        f.write(str(model)) # writes the changed model to new_model file



# #1. Bionetge run : average PZAP   over n run and obseravble is in the (p+1) column

    T,avg_ZAP,avg_SYK=run_bionetgen(n,p,p1,dir_name) #SYK
    T=np.asarray(T)
    avg_ZAP=np.asarray(avg_ZAP)#ZAP signa
    avg_SYK=np.asarray(avg_SYK)#SYK signal

    #2.Interpolate pSYKP  from 600 points to 2000 points  
    tnew,PZAP=spline(T, avg_ZAP, N,tstart) 
    tnew,PSYK=spline(T, avg_SYK, N,tstart)

    tnew=np.asarray(tnew)
    PZAP=np.asarray(PZAP) #total number of PZAP molecule in the simulation box of size Vc
    PSYK=np.asarray(PSYK) #total number of PZAP molecule in the simulation box of size Vc) 

    PZAP=PZAP/(Vc*z) #pZAP in uM unit
    PSYK=PSYK/(Vc*z)
    
    return tnew,PZAP,PSYK,N,dir_name



def calculate_residue(n,tnew,PZAP,PSYK,exp_time,exp_data,C1,C2,g,N,dir_name):
    
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
    ca_PZAP,h_PZAP=solve(tnew,N,y0z,PZAP,C1,C2,g)
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
    print(f"SSR={SSR}")

    
    return SSR    