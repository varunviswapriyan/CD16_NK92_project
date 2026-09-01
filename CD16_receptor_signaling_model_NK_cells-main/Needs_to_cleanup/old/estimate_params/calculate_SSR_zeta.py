import os, sys, shutil
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import bionetgen

from calculate_pZAP import run_bionetgen
from interpolation import spline





def zeta_PZAP(n,K1,K2,a0,ZAP0,SYK0,dir_name):


    p=1 # which column index which obseravable in bionetgen file for total ZAP70
    p1=2 # which column index which obseravable in bionetgen file for total SYK


    N=2000 #Ca signal time points and SSR at N timepoints

    tstart=0.0 # Fit to be start from which timepoint : interpolation of pZAP70 signal starts at 0 sec

    Vc=25.0 #pZAP molecules in the simulation box of size 25 um^3
    z=602.0 #constant factor to convert from molecules/um3 to uM
    
    print("k1 = "+str(K1))
    print("k2 = "+str(K2))
    print("A0= "+str(a0))
    print("ZAP0= "+str(ZAP0))
    print("SYK0= "+str(SYK0))


    dir_name = "zeta_runs/" + str(dir_name)
    model = bionetgen.bngmodel("JJ_zeta_0_July2.bngl")
    model.parameters.Kab = K1 # setting parameter kon to 1
    model.parameters.KU = K2 # setting parameter koff to 1
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



