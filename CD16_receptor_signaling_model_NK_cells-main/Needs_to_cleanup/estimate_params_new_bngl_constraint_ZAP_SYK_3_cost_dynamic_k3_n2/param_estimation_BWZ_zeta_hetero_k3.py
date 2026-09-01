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

def process_particle_zeta(args):
    """Helper function to process zeta PZAP calculation for a single particle."""
    i, param, outdir = args   
 
    a0_zeta = param[i, 3]
    ZAP0 = param[i, 4]
    SYK0 = param[i, 5]

    dir_name = str(outdir) + str(i)
    n = 3  # Number of bionetgen runs

    try:
        print(f'Processing zeta PZAP for particle {i}, directory_name={dir_name}')
        tnew_zeta, PZAP_zeta, PSYK_zeta, N, dir_name_zeta = zeta_PZAP(n, a0_zeta, ZAP0, SYK0, dir_name + "_zeta")
        return i, tnew_zeta, PZAP_zeta, PSYK_zeta, N, dir_name_zeta
    except Exception as e:
        print(f"Error processing zeta PZAP for particle {i}: {e}")
        return i, None, None, None, None, None

def process_particle_mixed(args):
    """Helper function to process mixed PZAP calculation for a single particle."""
    i, param, outdir = args

    a0_zeta = param[i, 3]
    a0_mixed = a0_zeta * 1.6
    ZAP0 = param[i, 4]
    SYK0 = param[i, 5]

    dir_name = str(outdir) + str(i)
    n = 3  # Number of bionetgen runs

    try:
        print(f'Processing mixed PZAP for particle {i}, directory_name={dir_name}')
        tnew_mixed, PZAP_mixed, PSYK_mixed, N, dir_name_mixed = mixed_PZAP(n, a0_mixed, ZAP0, SYK0, dir_name + "_mixed")
        return i, tnew_mixed, PZAP_mixed, PSYK_mixed, N, dir_name_mixed
    except Exception as e:
        print(f"Error processing mixed PZAP for particle {i}: {e}")
        return i, None, None, None, None, None
    

def process_particle_gamma(args):
    """Helper function to process gamma PZAP calculation for a single particle."""
    i, param, outdir = args

    a0_zeta = param[i, 3]
    a0_gamma = a0_zeta * 0.53
    ZAP0 = param[i, 4]
    SYK0 = param[i, 5]

    dir_name = str(outdir) + str(i)
    n = 3  # Number of bionetgen runs

    try:
        print(f'Processing gamma PZAP for particle {i}, directory_name={dir_name}')
        tnew_gamma, PZAP_gamma, PSYK_gamma, N, dir_name_gamma = gamma_PZAP(n, a0_gamma, ZAP0, SYK0, dir_name + "_gamma")
        return i, tnew_gamma, PZAP_gamma, PSYK_gamma, N, dir_name_gamma
    except Exception as e:
        print(f"Error processing gamma PZAP for particle {i}: {e}")
        return i, None, None, None, None, None



def ca_cost_func(param):
    param = 10 ** param 

    n_particles = param.shape[0]
    ssr_total = np.zeros(n_particles)
    
    #Number of particles and CPUs
    n_cpus = multiprocessing.cpu_count()
    print(f"number of CPUs available: {n_cpus}")
    
    # Prepare arguments for parallel processing
    args_zeta = [(i, param, outdir) for i in range(n_particles)]
    args_mixed = [(i, param, outdir) for i in range(n_particles)]
    args_gamma = [(i, param, outdir) for i in range(n_particles)]

    # Calculate PZAP signals in parallel using pools
    print("Calculating zeta PZAP in parallel...")
    with Pool(processes=multiprocessing.cpu_count()) as pool:
        zeta_results = pool.map(process_particle_zeta, args_zeta)
    
    print("Calculating mixed PZAP in parallel...")
    with Pool(processes=multiprocessing.cpu_count()) as pool:
        mixed_results = pool.map(process_particle_mixed, args_mixed)

    print("Calculating gamma PZAP in parallel...")
    with Pool(processes=multiprocessing.cpu_count()) as pool:
        gamma_results = pool.map(process_particle_gamma, args_gamma)

    # Calculate residues outside the pool
    print("Calculating residues...")
    for i in range(n_particles):
        try:
            C1 = param[i, 0]
            C2 = param[i, 1]
            g = param[i, 2]
            
            # Get zeta results
            zeta_i, tnew_zeta, PZAP_zeta, PSYK_zeta, N_zeta, dir_name_zeta = zeta_results[i]
            # Get mixed results  
            mixed_i, tnew_mixed, PZAP_mixed, PSYK_mixed, N_mixed, dir_name_mixed = mixed_results[i]
            # Get gamma results  
            gamma_i, tnew_gamma, PZAP_gamma, PSYK_gamma, N_gamma, dir_name_gamma = gamma_results[i]

            if tnew_zeta is None or tnew_mixed is None or tnew_gamma is None:
                # One of the PZAP calculations failed
                ssr_total[i] = float('inf')
                continue

            n=3    
            k3=(np.max(PZAP_zeta) + np.max(PZAP_mixed)) / 2.0  # Example calculation for k3, adjust as needed
            # Calculate residues outside the pool
            residue_zeta = calculate_residue(n, tnew_zeta, PZAP_zeta, PSYK_zeta, exp_time_zeta, exp_data_zeta, C1, C2, g, N_zeta, dir_name_zeta,k3)
            residue_mixed = calculate_residue(n, tnew_mixed, PZAP_mixed, PSYK_mixed, exp_time_mixed, exp_data_mixed, C1, C2, g, N_mixed, dir_name_mixed,k3)
            residue_gamma = calculate_residue(n, tnew_gamma, PZAP_gamma, PSYK_gamma, exp_time_gamma, exp_data_gamma, C1, C2, g, N_gamma, dir_name_gamma,k3)

            # Combine residues
            ssr_total[i] = residue_zeta + residue_mixed + residue_gamma
            
        except Exception as e:
            print(f"Error calculating residue for particle {i}: {e}")
            ssr_total[i] = float('inf')  # Assign a high cost for failed particles

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
    lb=[3.5,  -1.0,   -4.0,  +1.0, +2.5, +2.5] #lower bounds for C1,C2,g, A0_zeta, ZAP0, SYK
    ub=[6.0,   1.0,   -1.0,  +3.0,  3.2, 3.2] #upper bounds for C1,C2,g, A0_zeta, ZAP0, SYK

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
    
   

    
