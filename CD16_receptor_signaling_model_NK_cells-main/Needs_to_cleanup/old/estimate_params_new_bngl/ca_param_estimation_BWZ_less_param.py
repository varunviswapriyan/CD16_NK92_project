#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 21:33:25 2022

@author: ixn004
"""
#July 9, 2025 #ran the gamma code with single curve
#June 12, 2025
#April 5, 2025: trying all curves together using multiprocessing
#This code is used to estimate the parameters of the model using PSO: test for 1 curve
#April 1st 2025 - This is to test for 5 curves
#6th NOV 2023
#11th May 2023- run mean
#8th May, 2023- Mice run remove 3 parameters
#20th April, 2023 - Koff [0,1], bootstraap
#19th April. 2023 giving run for constant c1, c2,g
#14th September - thread lock is used 
#14th september - code is ready for parallel run
#All the functions checked on 5th September
# Ca ODE is checked on 6th September, where if I make PZAP signal is a unit vector, it matches with the the PZAP signal using MATLAB ode23
# here is the link to include link for CaODE in python https://apmonitor.com/pdc/index.php/Main/SolveDifferentialEquations
#7th September the objective function is considered after 25 sec


import os, sys, shutil
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import threading
import time
import multiprocessing

#from pyswarm import pso
from pyswarms.single.global_best import GlobalBestPSO #Date: Macrh 30,2024
from calculate_SSR_gamma import gamma_PZAP, calculate_residue
from calculate_SSR_zeta import zeta_PZAP    
from calculate_SSR_mixed import mixed_PZAP
from multiprocessing import Pool
import sys
import matplotlib.pyplot as plt

def process_particle_gamma(args):
    """Helper function to process a single particle."""
    i, param, outdir, exp_time_gamma, exp_data_gamma, exp_time_zeta, exp_data_zeta, exp_time_mixed, exp_data_mixed = args
    
    K1 = 0.0026675475345120147
    K2 = 0.1397364395199698



    C1 = param[i, 0]
    C2 = param[i, 1]
    g = param[i, 2]
    a0_zeta= param[i, 3]
    a0_mixed= param[i, 4]
    a0_gamma= param[i, 5]


    ZAP0 = 1131.0
    SYK0 = 754.0

    dir_name = str(outdir) + str(i)
    n = 3  # Number of bionetgen runs

    try:
        print(f'directory_name={dir_name}')
        tnew_gamma, PZAP_gamma, PSYK_gamma, N, dir_name_gamma = gamma_PZAP(n, K1, K2, a0_gamma, ZAP0, SYK0, dir_name)
        residue_gamma = calculate_residue(n, tnew_gamma, PZAP_gamma, PSYK_gamma, exp_time_gamma, exp_data_gamma, C1, C2, g, N, dir_name_gamma)

        total_residue = residue_gamma #+ residue_zeta + residue_mixed
    except Exception as e:
        print(f"Error processing particle {i}: {e}")
        total_residue = float('inf')  # Assign a high cost for failed particles

    return total_residue

def process_particle_zeta(args):
    """Helper function to process a single particle."""
    i, param, outdir, exp_time_gamma, exp_data_gamma, exp_time_zeta, exp_data_zeta, exp_time_mixed, exp_data_mixed = args   
    K1 = 0.0026675475345120147
    K2 = 0.1397364395199698



    C1 = param[i, 0]
    C2 = param[i, 1]
    g = param[i, 2]
    a0_zeta= param[i, 3]
    a0_mixed= param[i, 4]
    a0_gamma= param[i, 5]


    ZAP0 = 1131.0
    SYK0 = 754.0

    dir_name = str(outdir) + str(i)
    n = 3  # Number of bionetgen runs

    try:
        print(f'directory_name={dir_name}')
        tnew_zeta, PZAP_zeta, PSYK_zeta, N, dir_name_zeta = zeta_PZAP(n, K1, K2, a0_zeta, ZAP0, SYK0, dir_name)
        residue_zeta = calculate_residue(n, tnew_zeta, PZAP_zeta, PSYK_zeta, exp_time_zeta, exp_data_zeta, C1, C2, g, N, dir_name_zeta)

        total_residue = residue_zeta
    except Exception as e:
        print(f"Error processing particle {i}: {e}")
        total_residue = float('inf')  # Assign a high cost for failed particles

    return total_residue

def process_particle_mixed(args):
    """Helper function to process a single particle."""
    i, param, outdir, exp_time_gamma, exp_data_gamma, exp_time_zeta, exp_data_zeta, exp_time_mixed, exp_data_mixed = args
    K1 = 0.0026675475345120147
    K2 = 0.1397364395199698
    C1 = param[i, 0]
    C2 = param[i, 1]
    g = param[i, 2]
    a0_zeta= param[i, 3]
    a0_mixed= param[i, 4]
    a0_gamma= param[i, 5]

    ZAP0 = 1131.0
    SYK0 = 754.0


    dir_name = str(outdir) + str(i)
    n = 3  # Number of bionetgen runs

    try:
        print(f'directory_name={dir_name}')

        tnew_mixed, PZAP_mixed, PSYK_mixed, N, dir_name_mixed = mixed_PZAP(n, K1, K2, a0_mixed, ZAP0, SYK0, dir_name)
        residue_mixed = calculate_residue(n, tnew_mixed, PZAP_mixed, PSYK_mixed, exp_time_mixed, exp_data_mixed, C1, C2, g, N, dir_name_mixed)

        total_residue = residue_mixed
    except Exception as e:
        print(f"Error processing particle {i}: {e}")
        total_residue = float('inf')  # Assign a high cost for failed particles

    return total_residue

def ca_cost_func(param):
    param = 10 ** param 
    # print(param.shape)#dimension (particle_number, dimensions)
    # print(len(param))
    # print(f'param={param}')  
    # print(param.shape[0]) #number of particles

    args = [(i, param, outdir, exp_time_gamma, exp_data_gamma, exp_time_zeta, exp_data_zeta, exp_time_mixed, exp_data_mixed) for i in range(param.shape[0])]
    #print(f'args={args}') #list of tuples (i, param, outdir, exp_time_1, exp_data_1) and the length is equal to number of particles
  

    #Number of particles and CPUs
    n_cpus = multiprocessing.cpu_count()
    print(f"number of CPUs available: {n_cpus}")
    # Use multiprocessing to parallelize the loop
    with Pool(processes=multiprocessing.cpu_count()) as pool:
        ssr_gamma = pool.map(process_particle_gamma, args)

    with Pool(processes=multiprocessing.cpu_count()) as pool:
        ssr_zeta = pool.map(process_particle_zeta, args)

    with Pool(processes=multiprocessing.cpu_count()) as pool:
        ssr_mixed = pool.map(process_particle_mixed, args)

    ssr_total= np.array(ssr_gamma) + np.array(ssr_zeta) + np.array(ssr_mixed)

    return ssr_total



def main():
#5. PSO
    
    #1.start by importing your data  and bootstrapping  
    begin_time=time.time()
    global outdir 
    
    outdir='analysis' #this will contain the string analysis${i} output directories
    #os.mkdir(outdir)

    number_of_paricles=int(sys.argv[1]) #this will contain the string analysis${i} output directories
    iteration_number=int(sys.argv[2]) #this will contain the string analysis${i} output directories
    #2.start by importing mean data from Oscar

    data_gamma=pd.read_csv('../../../hCD16_Ca_expt/ca_data/scaled_data/gamma_scaled.dat',sep="\t", comment='#', header=None)
    data_zeta=pd.read_csv('../../../hCD16_Ca_expt/ca_data/scaled_data/zeta_scaled.dat',sep="\t", comment='#', header=None)
    data_mixed=pd.read_csv('../../../hCD16_Ca_expt/ca_data/scaled_data/mixed_scaled.dat',sep="\t", comment='#', header=None)

    global exp_time_gamma, exp_data_gamma, exp_time_zeta, exp_data_zeta, exp_time_mixed, exp_data_mixed
    #if the data has a time column:
    start_time=time.time()

    data_gamma=np.asarray(data_gamma) 
    data_zeta=np.asarray(data_zeta)
    data_mixed=np.asarray(data_mixed)
    
    exp_time_gamma=data_gamma[:,0]
    exp_data_gamma=data_gamma[:,1]

    exp_time_zeta=data_zeta[:,0]
    exp_data_zeta=data_zeta[:,1]

    exp_time_mixed=data_mixed[:,0]
    exp_data_mixed=data_mixed[:,1]

    #including adaptor: 6 parameters
    # parameter bounds [C1,C2,g, A0_zeta, A0_mixed, A0_gamma, ZAP0, SYK0]
    lb=[5.0,  -1.0,   -4.0,  +1.7, +1.7, +1.7] #lower bounds for C1,C2,g, A0_zeta, A0_mixed, A0_gamma
    ub=[6.0,   1.0,   -1.0,  3.0, 3.0, 3.0] #upper bounds for C1,C2,g, A0_zeta, A0_mixed, A0_gamma, ZAP0


    optimizer = GlobalBestPSO(n_particles=number_of_paricles, dimensions=6, options={'c1': 1.5, 'c2': 1.5, 'w': 0.5}, bounds=(lb, ub))
    #print(type(optimizer)

    # # Perform optimization
    residue, optimized_param = optimizer.optimize(ca_cost_func, iters=iteration_number)
    # #define a lock aloow one thread at a time, all other thread muct wait until the lock is released

    print(f'shape of residue {residue.shape}')
    print(f'shape of optimized_param {optimized_param.shape}')

  

    thlock=threading.Lock()       
    
    print(optimized_param)
    print(residue)
    residue_small=residue**(1/4)
    print(f' residue ** (1/4) = {residue_small}')
    xx=math.sqrt(residue)
    yy=math.sqrt(xx) 
    
    print(type(optimized_param))
    optimized_param=np.asarray(optimized_param)
    
    f = open(str(outdir)+"_param_residue.dat", 'w')
    f.write("The optimized parameter is "+str(optimized_param)+"\n Residue = "+str(residue)+"\n Squared residue = " +str(xx)+"\n 1/4th residue = "+str(yy))
    f.close()
        
    thlock.acquire() #lock on
    f2 = open("common.dat", 'a') #append in a common file
    f2.write(str(outdir)+"\t"+str(optimized_param)+"\n")
    thlock.release() #lock off
    

                
    thlock.acquire() #lock on
    f3 = open("common_parameters.dat", 'a')
    f3.write(str(optimized_param[0])+"\t"+str(optimized_param[1])+"\t"+str(optimized_param[2])+"\t"+str(optimized_param[3])+"\t"+str(optimized_param[4])+"\t"+str(optimized_param[5])+"\t"+str(yy)+"\n")
    thlock.release() #lock off

    end_time=time.time()
    time_taken=(end_time-start_time)/60.0
    print(f'The time taken for the process is {time_taken} min')
    print(f'The process ends now')
        

if __name__=="__main__" :
    main()    
    
   

    
