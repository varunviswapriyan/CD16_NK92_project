import os, sys, shutil, functools, hashlib, glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import scipy.stats
import bionetgen
import threading

#from pyswarm import pso
from pyswarms.single.global_best import GlobalBestPSO
from run_bngl import run_bionetgen
from interplation import spline
from ca_ODE import calcium, solve
from interpolation_exp_data import exp_idata
from oscar_data_import import expt_import
from ca_bootstrap import bootstrap_ca
from data_slice import data_after_tstart
import sys


    
def ZAP_only(n,K1,K2,C1,C2,g,exp_time,exp_data,dir_name):
    
    
            
            p=1 # which column index which obseravable in bionetgen file
            
            N=2000 #Ca signal time points and SSR at N timepoints
            
            tstart=26.0 # Fit to be start from which timepoint : interpolation of pZAP70 signal starts at 25 sec
            
            Vc=25.0 #pZAP molecules in the simulation box of size 25 um^3
            z=602.0 #constant factor to convert from molecules/um3 to uM
            
            # print("k1 = "+str(K1))
            # print("k2 = "+str(K2))
            # print("C1 = "+str(C1))
            # print("C2 = "+str(C2))
            # print("g = "+str(g))
        

        #0. Bionetgen parameter set 
            model = bionetgen.bngmodel("zeta_0.bngl")
            model.parameters.Kab = K1 # setting parameter kon to 1
            model.parameters.KU = K2 # setting parameter koff to 1
            
            #print(model)

            #print model in directiry name_gamma_HPC.bngl    5 folders
            with open(dir_name+"_zeta.bngl", "w") as f:
                f.write(str(model)) # writes the changed model to new_model file


        
        # #1. Bionetge run : average PZAP   over n run and obseravble is in the (p+1) column
        
            T,avg=run_bionetgen(n,p,dir_name) 
            T=np.asarray(T)
            avg=np.asarray(avg)
            

            # print(T)
            # print(avg)
            # print(len(T))
            
            
        #2.Interpolate PZAP  from 600 points to 2000 points   
            tnew,PZAP=spline(T, avg, N,tstart)
            tnew=np.asarray(tnew)
            PZAP=np.asarray(PZAP) #total number of PZAP molecule in the simulation box of size Vc
            PZAP=PZAP/(Vc*z) #pZAP in uM unit
            # plt.plot(T,avg,'g*',tnew,PZAP,'r-') 
            # plt.show()  
            
        
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

            
        #5. Ca Signal   from the ODE model uisng the PZAP signal  
            #CA0=33212 #do not change 
            y0 = [CA0, 1]
            ca,h=solve(tnew,N,y0,PZAP,C1,C2,g)
            ca=np.asarray(ca)
            h=np.asarray(h)
            # plt.plot(tnew,ca,'m')
            # plt.show()
      
            
                
        #6. fit starts from tstart, so, MODEL and EXPT data containts data after tstart    
            MODEL=[]
            for i in range (count,len(tnew),1): 
                MODEL.append(ca[i])
        
            TT=np.asarray(TT)
            EXPT=np.asarray(EXPT)
            MODEL=np.asarray(MODEL)
            #print(TT)
            #print(TT)
            
        #7. Finding objective function only at t>25 sec
            ca_diff=(idata-ca)
            OBJ=(EXPT-MODEL)
            
            #print(ca_diff)
            SSR=np.sum(np.square(OBJ))
            print(SSR)


            
            # plot_data=np.column_stack([tnew,ca,idata])
            # file = open('plot_fit.dat', 'w')
            # np.savetxt("plot_fit.dat" , plot_data, fmt=['%f','%f','%f'])
            #plt.plot(TT,MODEL,'g-',tnew,idata,'r-',TT,OBJ,'b--',time,exp_data,'m')
            #plt.show()
            
            
     
            
            return SSR,ca,tnew   