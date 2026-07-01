"""
Stretching experiment of LiH - Parallelised Version
"""

import os
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count
from qchem.utils import (
    MoleculeData,
    DMDMWorkflow,
)
from dmdm.interface import DMDM
import qrunch as qc

def run_single_factor(args):
    """Worker function for each factor - must be at module level for pickling"""
    factor, num_active_orbitals, num_active_electrons, num_states = args
    
    print(f"Processing factor: {factor}")
    
    # Create calculator fresh in each worker (avoid pickling issues)
    calculator = (qc.calculator_creator()
                    .vqe()
                    .iterative()
                    .standard()
                    .with_options(
                            options=qc.options.IterativeVqeOptions(max_iterations=500)
                    )
                    .create())
    
    mol = [('Li', 0.0, 0.0, 0.0), ('H', 0.0, 0.0, factor)]
    
    workflow = DMDMWorkflow(
        basis="cc-pvdz",
        molecule=mol,
        num_active_orbitals=num_active_orbitals,
        num_active_electrons=num_active_electrons,
        num_states=num_states,
        
        calculator=calculator,
        casci_like=True
    )
    
    # Run all calculations
    result_vqe = workflow.run_quantum_vqe()
    result_casci = workflow.run_classical_casci()
    
    return {
        'factor': factor,
        'vqe': {
            "state": [*range(1, len(result_vqe["exc_energies_ev"])+1)],
            "energy_ev": list(result_vqe["exc_energies_ev"]),
            "oscillator_strength": list(result_vqe["oscillator_strengths"])
        },
        'casci': {
            "state": [*range(1, len(result_casci["exc_energies_ev"])+1)],
            "energy_ev": list(result_casci["exc_energies_ev"]),
            "oscillator_strength": list(result_casci["oscillator_strengths"])
        },
    }


if __name__ == "__main__":
    # Define parameters
    factors = (
        np.round(np.linspace(0.5, 1.5, 11), 2).tolist() + 
        [1.595] + 
        np.round(np.linspace(1.6, 4.0, 25), 2).tolist()
    )
    num_workers = min(cpu_count(), 24)  # Limit workers to avoid memory issues
    
    # Prepare arguments for parallel execution
    args_list = [
        (factor, 6, 4, 40)  # num_active_orbitals, num_active_electrons, num_states
        for factor in factors
    ]
    
    # Execute in parallel
    with Pool(processes=num_workers) as pool:
        results = pool.map(run_single_factor, args_list)
    
    # Consolidate results (same structure as original)
    vqe_dfs = []
    casci_dfs = []
    
    for result in results:
        factor = result['factor']
        
        vqe_dfs.append(pd.DataFrame({
            "state": result['vqe']['state'],
            "energy_ev": result['vqe']['energy_ev'],
            "stretch": [factor] * len(result['vqe']['state']),
            "oscillator_strength": result['vqe']['oscillator_strength']
        }))
        
        casci_dfs.append(pd.DataFrame({
            "state": result['casci']['state'],
            "energy_ev": result['casci']['energy_ev'],
            "stretch": [factor] * len(result['casci']['state']),
            "oscillator_strength": result['casci']['oscillator_strength']
        }))
        
    
    # Save results
    os.makedirs("lih_stretching", exist_ok=True)
    
    vqe_df = pd.concat(vqe_dfs, axis=0).reset_index(drop=True)
    casci_df = pd.concat(casci_dfs, axis=0).reset_index(drop=True)
    
    vqe_df.to_csv("lih_stretching/vqe_lih_stretching.csv")
    casci_df.to_csv("lih_stretching/casci_lih_stretching.csv")
    
    print(f"Parallel processing complete. Results saved to lih_stretching/")