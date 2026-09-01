#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on August 25, 2025

Parameter estimation for 5 parameters: C1, C2, g, k3, k4
Cost function: residue_zeta + residue_mixed + residue_gamma

@author: ixn004
"""

import os, sys, shutil
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import threading
import time
import multiprocessing

from pyswarms.single.global_best import GlobalBestPSO
from ZAP_SYK_SSR_fit import calculate_residue


def ca_cost_func(param):
    """
    Cost function for 4-parameter optimization: C1, C2, g, k4
    k3 is fixed at 0.06
    Uses pre-loaded PZAP/PSYK data from files instead of calculating them.

    """
    
    # Convert from log10 space for 4 parameters
    param_converted = 10 ** param  # Convert C1, C2, g, k4 from log space

    n_particles = param.shape[0]
    ssr_total = np.zeros(n_particles)
    
    print(f"Calculating residues for {n_particles} particles...")
    
    # Calculate residues for each particle
    for i in range(n_particles):
        try:
            C1 = param_converted[i, 0]
            C2 = param_converted[i, 1]
            g = param_converted[i, 2]
            k3 = 0.08  # Fixed value for k3
            k4 = param_converted[i, 3]
            
            n = 3  # Hill coefficient
            N = 2000  # Number of time points for calcium ODE
            
            # Calculate residues for each mechanism using pre-loaded data
            # calculate_residue returns (SSR, TT, EXPT, MODEL) - we only need SSR
            residue_zeta, _, _, _ = calculate_residue(n, tnew_zeta, pZAP_zeta, pSYK_zeta, 
                                                    exp_time_zeta, exp_data_zeta, 
                                                    C1, C2, g, N, k3, k4)
            
            residue_mixed, _, _, _ = calculate_residue(n, tnew_mixed, pZAP_mixed, pSYK_mixed, 
                                                     exp_time_mixed, exp_data_mixed, 
                                                     C1, C2, g, N, k3, k4)
            
            residue_gamma, _, _, _ = calculate_residue(n, tnew_gamma, pZAP_gamma, pSYK_gamma, 
                                                     exp_time_gamma, exp_data_gamma, 
                                                     C1, C2, g, N, k3, k4)

            # Combined cost function: residue_zeta + residue_mixed + residue_gamma
            ssr_total[i] = residue_zeta + residue_mixed + residue_gamma
            
            if i % 5 == 0:  # Print every 5th particle to reduce output
                print(f"Particle {i}: C1={C1:.2e}, C2={C2:.2e}, g={g:.2e}, k3={k3:.2e} (fixed), k4={k4:.2e}")
                print(f"  Residues: zeta={residue_zeta:.4f}, mixed={residue_mixed:.4f}, gamma={residue_gamma:.4f}")
                print(f"  Total cost: {ssr_total[i]:.4f}")
            
        except Exception as e:
            print(f"Error calculating residue for particle {i}: {e}")
            ssr_total[i] = float('inf')  # Assign a high cost for failed particles

    return ssr_total


def main():
    """Main function for parameter estimation using PSO."""
    
    begin_time = time.time()
    
    # Get command line arguments
    if len(sys.argv) < 3:
        print("Usage: python estimate_k3_k4.py <number_of_particles> <number_of_iterations>")
        print("Example: python estimate_k3_k4.py 20 50")
        sys.exit(1)
        
    number_of_particles = int(sys.argv[1])
    iteration_number = int(sys.argv[2])
    
    print(f"Starting parameter estimation with {number_of_particles} particles and {iteration_number} iterations")
    print("k3 is fixed at 0.06 - optimizing C1, C2, g, k4")

    # Load experimental data
    print("Loading experimental data...")
    data_gamma = pd.read_csv('../../../hCD16_Ca_expt/ca_data/scaled_data/gamma_scaled.dat', 
                            sep="\t", comment='#', header=None)
    data_zeta = pd.read_csv('../../../hCD16_Ca_expt/ca_data/scaled_data/zeta_scaled.dat', 
                           sep="\t", comment='#', header=None)
    data_mixed = pd.read_csv('../../../hCD16_Ca_expt/ca_data/scaled_data/mixed_scaled.dat', 
                            sep="\t", comment='#', header=None)

    # Make experimental data global for access in cost function
    global exp_time_gamma, exp_data_gamma, exp_time_zeta, exp_data_zeta, exp_time_mixed, exp_data_mixed
    
    data_gamma = np.asarray(data_gamma) 
    data_zeta = np.asarray(data_zeta)
    data_mixed = np.asarray(data_mixed)
    
    exp_time_gamma = data_gamma[:, 0]
    exp_data_gamma = data_gamma[:, 1]
    exp_time_zeta = data_zeta[:, 0]
    exp_data_zeta = data_zeta[:, 1]
    exp_time_mixed = data_mixed[:, 0]
    exp_data_mixed = data_mixed[:, 1]

    print("Experimental data loaded successfully")
    
    # Load pZAP and pSYK data from files
    print("Loading pZAP and pSYK simulation data...")
    
    # Make simulation data global for access in cost function
    global tnew_zeta, pZAP_zeta, pSYK_zeta, tnew_mixed, pZAP_mixed, pSYK_mixed, tnew_gamma, pZAP_gamma, pSYK_gamma
    
    # Load zeta data
    zeta_data = pd.read_csv('zeta_results.txt', sep='\t', header=None)
    tnew_zeta = zeta_data.iloc[:, 0].values
    pZAP_zeta = zeta_data.iloc[:, 1].values
    pSYK_zeta = zeta_data.iloc[:, 2].values
    
    # Load mixed data
    mixed_data = pd.read_csv('mixed_results.txt', sep='\t', header=None)
    tnew_mixed = mixed_data.iloc[:, 0].values
    pZAP_mixed = mixed_data.iloc[:, 1].values
    pSYK_mixed = mixed_data.iloc[:, 2].values
    
    # Load gamma data
    gamma_data = pd.read_csv('gamma_results.txt', sep='\t', header=None)
    tnew_gamma = gamma_data.iloc[:, 0].values
    pZAP_gamma = gamma_data.iloc[:, 1].values
    pSYK_gamma = gamma_data.iloc[:, 2].values
    
    print(f"Simulation data loaded successfully:")
    print(f"  Zeta: {len(tnew_zeta)} time points")
    print(f"  Mixed: {len(tnew_mixed)} time points") 
    print(f"  Gamma: {len(tnew_gamma)} time points")

    # Parameter bounds for optimization (in log10 space for 4 parameters)
    # [log10(C1), log10(C2), log10(g), log10(k4)]
    # k3 is fixed at 0.06
    
    # Bounds for optimization parameters (log10 space)
    lb = [2.0,   # log10(C1) - lower bound around 1e2
          -2.0,  # log10(C2) - lower bound around 1e-2  
          -4.0,  # log10(g) - lower bound around 1e-4
          -3.0]  # log10(k4) - lower bound around 1e-4

    ub = [6.0,   # log10(C1) - upper bound around 1e6
          1.0,   # log10(C2) - upper bound around 10
          -1.0,  # log10(g) - upper bound around 0.1
          +1.0]   # log10(k4) - upper bound around 0.16

    print(f"Parameter bounds (4 parameters in log10 space, k3 fixed at 0.06):")
    print(f"Lower bounds: {lb}")
    print(f"Upper bounds: {ub}")

    # Initialize PSO optimizer
    optimizer = GlobalBestPSO(n_particles=number_of_particles, 
                              dimensions=4,  # Only 4 parameters to optimize (C1, C2, g, k4)
                              options={'c1': 1.5, 'c2': 1.5, 'w': 0.5}, 
                              bounds=(lb, ub))

    print("Starting PSO optimization...")
    start_optimization = time.time()

    # Perform optimization
    residue, optimized_param = optimizer.optimize(ca_cost_func, iters=iteration_number)

    end_optimization = time.time()
    optimization_time = (end_optimization - start_optimization) / 60.0
    
    print(f"\nOptimization completed in {optimization_time:.2f} minutes")
    print(f'Residue shape: {residue.shape}')
    print(f'Optimized parameters shape: {optimized_param.shape}')

    # Extract optimized parameters (convert from log10 space)
    C1_opt = 10**optimized_param[0]
    C2_opt = 10**optimized_param[1]
    g_opt = 10**optimized_param[2]
    k3_opt = 0.08  # Fixed value
    k4_opt = 10**optimized_param[3]
    
    print("\n" + "="*60)
    print("OPTIMIZATION RESULTS:")
    print("="*60)
    print(f"Optimal parameters:")
    print(f"C1_optimal = {C1_opt:.4e}")
    print(f"C2_optimal = {C2_opt:.4e}")  
    print(f"g_optimal = {g_opt:.4e}")
    print(f"k3_optimal = {k3_opt:.4e} (fixed)")
    print(f"k4_optimal = {k4_opt:.4e}")
    print(f"Final residue: {residue:.6f}")
    
    residue_small = residue**(1/4)
    print(f'Residue^(1/4) = {residue_small:.6f}')
    
    # Thread lock for file operations
    thlock = threading.Lock()       
    
    # Save results to files
    thlock.acquire()
    try:
        # Save detailed results
        with open(f"analysis_param_residue.dat", 'w') as f:
            f.write(f"Optimization Results - 4 Parameters (C1, C2, g, k4) with k3 fixed at 0.06\n")
            f.write(f"Number of particles: {number_of_particles}\n")
            f.write(f"Number of iterations: {iteration_number}\n")
            f.write(f"Optimization time: {optimization_time:.2f} minutes\n\n")
            f.write(f"Optimized parameters (log10 space): {optimized_param}\n")
            f.write(f"C1_optimal = {C1_opt:.6e}\n")
            f.write(f"C2_optimal = {C2_opt:.6e}\n")
            f.write(f"g_optimal = {g_opt:.6e}\n")
            f.write(f"k3_optimal = {k3_opt:.6e} (fixed)\n")
            f.write(f"k4_optimal = {k4_opt:.6e}\n")
            f.write(f"Final residue = {residue:.6f}\n")
            f.write(f"Residue^(1/4) = {residue_small:.6f}\n")
            
        # Append to common files
        with open("common.dat", 'a') as f2:
            f2.write(f"analysis\t{optimized_param}\n")
            
        with open("common_parameters.dat", 'a') as f3:
            f3.write(f"{optimized_param[0]:.6f}\t{optimized_param[1]:.6f}\t{optimized_param[2]:.6f}\t")
            f3.write(f"{np.log10(k3_opt):.6f}\t{optimized_param[3]:.6f}\t{residue_small:.6f}\n")
            
    finally:
        thlock.release()

    # Print final summary
    end_time = time.time()
    total_time = (end_time - begin_time) / 60.0
    
    print("\n" + "="*60)
    print("EXECUTION SUMMARY:")
    print("="*60)
    print(f'Total execution time: {total_time:.2f} minutes')
    print(f'Optimization time: {optimization_time:.2f} minutes')
    print(f'Setup and I/O time: {total_time - optimization_time:.2f} minutes')
    print(f'Cost function: residue_zeta + residue_mixed + residue_gamma')
    print(f'Final combined residue: {residue:.6f}')
    print(f'Results saved to: analysis_param_residue.dat')
    print('Process completed successfully!')
    print("="*60)


if __name__ == "__main__":
    main()
