#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ligand Concentration Variation Analysis for CD16 Receptor Signaling

This script varies ligand0 concentration from 60 to 500 and analyzes 
the maximum PZAP values for different mechanisms (zeta, mixed, gamma).

Created on August 29, 2025
@author: ixn004
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
import os
import sys

# Import required modules
from pyswarms.single.global_best import GlobalBestPSO
from calculate_pZAP import run_bionetgen
from interpolation import spline
from ca_ODE import calcium, solve
from interpolation_exp_data import exp_idata
from data_slice import data_after_tstart
from ZAP_SYK_SSR_fit import calculate_residue
from calculate_SSR_gamma import gamma_PZAP
from calculate_SSR_zeta import zeta_PZAP    
from calculate_SSR_mixed import mixed_PZAP


def setup_parameters():
    """Set up the simulation parameters."""
    # Model parameters
    n = 3
    
    # Receptor concentrations
    a0_zeta = 250.0
    a0_mixed = a0_zeta * 1.6 
    a0_gamma = a0_zeta * 0.53
    ZAP0 = 1131.0 
    SYK0 = 754.0
    
    # Optimized parameters
    C1_optimal = 3.0386e+05
    C2_optimal = 3.1280e+00
    g_optimal = 8.0032e-03
    k3_optimal = 6.6179e-02
    k4_optimal = 3.1675e-01
    
    # Other parameters
    dir_name = "analysis"
    
    return {
        'n': n,
        'a0_zeta': a0_zeta,
        'a0_mixed': a0_mixed,
        'a0_gamma': a0_gamma,
        'ZAP0': ZAP0,
        'SYK0': SYK0,
        'C1': C1_optimal,
        'C2': C2_optimal,
        'g': g_optimal,
        'k3': k3_optimal,
        'k4': k4_optimal,
        'dir_name': dir_name
    }


def print_parameters(params):
    """Print the simulation parameters."""
    print("SIMULATION PARAMETERS")
    print("=" * 50)
    print(f"C1 = {params['C1']:.4e}")
    print(f"C2 = {params['C2']:.4e}")
    print(f"g = {params['g']:.4e}")
    print(f"k3 = {params['k3']:.4e}")
    print(f"k4 = {params['k4']:.4e}")
    print()
    print(f"a0_zeta = {params['a0_zeta']}")
    print(f"a0_mixed = {params['a0_mixed']}")
    print(f"a0_gamma = {params['a0_gamma']}")
    print(f"ZAP0 = {params['ZAP0']}")
    print(f"SYK0 = {params['SYK0']}")
    print("=" * 50)


def run_ligand_variation_analysis(params, ligand_concentrations=None):
    """
    Run the ligand concentration variation analysis.
    
    Parameters:
    -----------
    params : dict
        Dictionary containing simulation parameters
    ligand_concentrations : list, optional
        List of ligand concentrations to test. If None, uses default range.
    
    Returns:
    --------
    dict : Dictionary containing max PZAP results
    """
    # Define default ligand concentrations if not provided
    if ligand_concentrations is None:
        ligand_concentrations = [60, 120, 180, 240, 300, 360, 420, 480, 540, 600]
    # Storage for maximum PZAP values
    max_PZAP_results = {
        'ligand0': [],
        'max_PZAP_zeta': [],
        'max_PZAP_mixed': [],
        'max_PZAP_gamma': []
    }
    
    print("\nStarting ligand concentration variation analysis...")
    print("=" * 60)
    
    for ligand0 in ligand_concentrations:
        print(f"\nProcessing ligand0 = {ligand0}")
        
        # Calculate PZAP and PSYK for each mechanism
        tnew_zeta, PZAP_zeta, PSYK_zeta, N, dir_name_zeta = zeta_PZAP(
            params['n'], params['a0_zeta'], params['ZAP0'], params['SYK0'], 
            ligand0, params['dir_name']
        )
        tnew_mixed, PZAP_mixed, PSYK_mixed, N, dir_name_mixed = mixed_PZAP(
            params['n'], params['a0_mixed'], params['ZAP0'], params['SYK0'], 
            ligand0, params['dir_name']
        )
        tnew_gamma, PZAP_gamma, PSYK_gamma, N, dir_name_gamma = gamma_PZAP(
            params['n'], params['a0_gamma'], params['ZAP0'], params['SYK0'], 
            ligand0, params['dir_name']
        )
        
        # Save results with ligand0 subscript
        save_results_to_files(ligand0, tnew_zeta, PZAP_zeta, PSYK_zeta,
                             tnew_mixed, PZAP_mixed, PSYK_mixed,
                             tnew_gamma, PZAP_gamma, PSYK_gamma)
        
        # Calculate maximum PZAP values
        max_PZAP_zeta = np.max(PZAP_zeta)
        max_PZAP_mixed = np.max(PZAP_mixed)
        max_PZAP_gamma = np.max(PZAP_gamma)
        
        # Store results
        max_PZAP_results['ligand0'].append(ligand0)
        max_PZAP_results['max_PZAP_zeta'].append(max_PZAP_zeta)
        max_PZAP_results['max_PZAP_mixed'].append(max_PZAP_mixed)
        max_PZAP_results['max_PZAP_gamma'].append(max_PZAP_gamma)
        
        print(f"  Max PZAP - Zeta: {max_PZAP_zeta:.2f}, Mixed: {max_PZAP_mixed:.2f}, Gamma: {max_PZAP_gamma:.2f}")
    
    print("\n" + "=" * 60)
    print("Ligand concentration variation analysis completed!")
    print(f"Processed {len(ligand_concentrations)} different concentrations")
    
    return max_PZAP_results


def save_results_to_files(ligand0, tnew_zeta, PZAP_zeta, PSYK_zeta,
                         tnew_mixed, PZAP_mixed, PSYK_mixed,
                         tnew_gamma, PZAP_gamma, PSYK_gamma):
    """Save simulation results to files with ligand0 subscript."""
    
    # Save zeta results
    zeta_filename = f"zeta_results_{ligand0}.txt"
    with open(zeta_filename, "w") as f:
        for i in range(len(tnew_zeta)):
            f.write(f"{tnew_zeta[i]}\t{PZAP_zeta[i]}\t{PSYK_zeta[i]}\n")
    
    # Save mixed results
    mixed_filename = f"mixed_results_{ligand0}.txt"
    with open(mixed_filename, "w") as f:
        for i in range(len(tnew_mixed)):
            f.write(f"{tnew_mixed[i]}\t{PZAP_mixed[i]}\t{PSYK_mixed[i]}\n")
    
    # Save gamma results
    gamma_filename = f"gamma_results_{ligand0}.txt"
    with open(gamma_filename, "w") as f:
        for i in range(len(tnew_gamma)):
            f.write(f"{tnew_gamma[i]}\t{PZAP_gamma[i]}\t{PSYK_gamma[i]}\n")
    
    print(f"  Saved: {zeta_filename}, {mixed_filename}, {gamma_filename}")


def create_visualizations(max_PZAP_results):
    """Create bar plots and trend lines for the analysis results."""
    
    # Convert to numpy arrays for easier manipulation
    ligand_conc = np.array(max_PZAP_results['ligand0'])
    max_zeta = np.array(max_PZAP_results['max_PZAP_zeta'])
    max_mixed = np.array(max_PZAP_results['max_PZAP_mixed'])
    max_gamma = np.array(max_PZAP_results['max_PZAP_gamma'])
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Bar plot 1: Individual mechanisms
    x_pos = np.arange(len(ligand_conc))
    width = 0.25
    
    bars1 = ax1.bar(x_pos - width, max_zeta, width, label='Zeta', alpha=0.8, color='blue')
    bars2 = ax1.bar(x_pos, max_mixed, width, label='Mixed', alpha=0.8, color='green')
    bars3 = ax1.bar(x_pos + width, max_gamma, width, label='Gamma', alpha=0.8, color='red')
    
    ax1.set_xlabel('Ligand Concentration')
    ax1.set_ylabel('Max PZAP')
    ax1.set_title('Maximum PZAP vs Ligand Concentration\n(Individual Mechanisms)')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(ligand_conc, rotation=45)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            ax1.annotate(f'{height:.0f}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=8)
    
    # Line plot 2: Trends
    ax2.plot(ligand_conc, max_zeta, 'o-', label='Zeta', linewidth=2, markersize=6, color='blue')
    ax2.plot(ligand_conc, max_mixed, 's-', label='Mixed', linewidth=2, markersize=6, color='green')
    ax2.plot(ligand_conc, max_gamma, '^-', label='Gamma', linewidth=2, markersize=6, color='red')
    
    ax2.set_xlabel('Ligand Concentration')
    ax2.set_ylabel('Max PZAP')
    ax2.set_title('Maximum PZAP vs Ligand Concentration\n(Trend Lines)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('max_PZAP_vs_ligand_concentration.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return ligand_conc, max_zeta, max_mixed, max_gamma


def print_summary_statistics(ligand_conc, max_zeta, max_mixed, max_gamma):
    """Print summary statistics of the analysis."""
    
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)
    print(f"{'Ligand':<8} {'Max_Zeta':<10} {'Max_Mixed':<11} {'Max_Gamma':<11}")
    print("-" * 45)
    for i in range(len(ligand_conc)):
        print(f"{ligand_conc[i]:<8} {max_zeta[i]:<10.1f} {max_mixed[i]:<11.1f} {max_gamma[i]:<11.1f}")
    
    print(f"\nRanges:")
    print(f"Zeta:  {np.min(max_zeta):.1f} - {np.max(max_zeta):.1f} (fold change: {np.max(max_zeta)/np.min(max_zeta):.2f})")
    print(f"Mixed: {np.min(max_mixed):.1f} - {np.max(max_mixed):.1f} (fold change: {np.max(max_mixed)/np.min(max_mixed):.2f})")
    print(f"Gamma: {np.min(max_gamma):.1f} - {np.max(max_gamma):.1f} (fold change: {np.max(max_gamma)/np.min(max_gamma):.2f})")


def save_results_to_csv(max_PZAP_results):
    """Save results to CSV file."""
    results_df = pd.DataFrame(max_PZAP_results)
    results_df.to_csv('max_PZAP_vs_ligand_concentration.csv', index=False)
    print(f"\nResults saved to: max_PZAP_vs_ligand_concentration.csv")
    print(f"Plot saved to: max_PZAP_vs_ligand_concentration.png")


def main():
    """Main function to run the complete analysis."""
    
    print("CD16 Receptor Signaling - Ligand Concentration Variation Analysis")
    print("=" * 70)
    
    # Setup parameters
    params = setup_parameters()
    print_parameters(params)
    
    # Option to use reduced set for testing (uncomment the line below)
    # ligand_concentrations = [60, 300, 500]  # For quick testing
    ligand_concentrations = None  # Use full range
    
    # Run analysis
    max_PZAP_results = run_ligand_variation_analysis(params, ligand_concentrations)
    
    # Create visualizations
    ligand_conc, max_zeta, max_mixed, max_gamma = create_visualizations(max_PZAP_results)
    
    # Print summary statistics
    print_summary_statistics(ligand_conc, max_zeta, max_mixed, max_gamma)
    
    # Save results
    save_results_to_csv(max_PZAP_results)
    
    print("\n" + "=" * 70)
    print("Analysis completed successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
