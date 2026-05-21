"""
Runs an experiment to compare VQE with CACI, CASSCF
"""

import argparse
import subprocess
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyscf import gto, scf, mcscf, ao2mo
from qchem.utils import  (
    MoleculeData,
    get_hf_gse_from_mol,
    one_electron_integral_transform,
    two_electron_integral_transform,
    DMDMWorkflow,
    CalculationMode,
    scale_molecule_v2,
    write_dalton_molecule_file,
    parse_dalton_output,
    gaussian,
)

from dmdm.interface import DMDM
import qrunch as qc

def main():


    parser = argparse.ArgumentParser()

    parser.add_argument(
        "molecule",
        type=str,
        help="molecule name"
    )

    parser.add_argument(
        "num_active_orbitals",
        type=int,
        help="num_active_orbitals"
    )

    parser.add_argument(
        "num_active_electrons",
        type=int,
        help="num_active_electrons"
    )

    parser.add_argument(
        "num_states",
        type=int,
        help="num_states"
    )

    parser.add_argument(
        "-o",
        type=str,
        help="Output folder path for experiment"
    )

    parser.add_argument(
        "-b",
        type=str,
        default="cc-pvdz",
        help="Molecule basis"
    )

    args = parser.parse_args()

    print("Running VQE vs. CASCI comparison\n")
    print(args)
    print("\n")

    # Create output folder
    Path(f"{args.o}").mkdir(parents=True, exist_ok=True)

    # Change the global working directory
    os.chdir(args.o)



    # Get molecule
    molecule = MoleculeData.molecules[args.molecule]["coords"]

    calculator = (
        qc.calculator_creator()
        .vqe()
        .iterative()
        .standard()
        .with_options(options=qc.options.IterativeVqeOptions(max_iterations=500))
        .choose_stopping_criterion()
        .patience(patience=10, threshold=1e-10)
        .create()
    )

    workflow = DMDMWorkflow(
        basis=args.b,
        molecule=molecule,
        num_active_orbitals=args.num_active_orbitals,
        num_active_electrons=args.num_active_electrons,
        num_states=args.num_states,
        mode=CalculationMode.BOTH,
        vqe_max_iterations=500,
        casci_like=True,
        calculator=calculator
    )

    # take the time
    start = time.time()
    workflow.run_quantum_vqe()
    vqe_end = time.time() - start
    print(f"VQE (cas-like) time (seconds): {vqe_end}")


    start = time.time()
    workflow.run_classical_casci()
    vqe_end = time.time() - start
    print(f"CASCI time (seconds): {vqe_end}")


    start = time.time()
    workflow.run_classical_casci_dmdm()
    vqe_end = time.time() - start
    print(f"CASCI time (dmdm) (seconds): {vqe_end}")


    calculator = (
        qc.calculator_creator()
        .vqe()
        .iterative_with_orbital_optimization()
        .standard()
        .with_options(options=qc.options.IterativeVqeOptions(max_iterations=500))
        .choose_stopping_criterion()
        .patience(patience=10, threshold=1e-10)
        .create()
    )


    workflow2 = DMDMWorkflow(
        basis=args.b,
        molecule=molecule,
        num_active_orbitals=args.num_active_orbitals,
        num_active_electrons=args.num_active_electrons,
        num_states=args.num_states,
        mode=CalculationMode.BOTH,
        vqe_max_iterations=500,
        casci_like=False,
        calculator=calculator
    )

    start = time.time()
    workflow2.run_quantum_vqe()
    vqe_end = time.time() - start
    print(f"VQE time (seconds): {vqe_end}")

        # Save the results as dfs
    pd.DataFrame(
        workflow._vqe_results
    ).to_csv(
        "vqe_casci_like.csv"
    )

    pd.DataFrame(
        workflow._casci_results
    ).to_csv(
        "casci.csv"
    )

    pd.DataFrame(
        workflow._casci_dmdm_results
    ).to_csv(
        "casci_dmdm.csv"
    )

    pd.DataFrame(
        workflow2._vqe_results
    ).to_csv(
        "vqe.csv"
    )

    # save the compute times
    pd.DataFrame(
        {
            "method": ["vqe_as", "vqe", "casci", "casci_dmdm"],
            "time": [
                workflow.vqe_time,
                workflow2.vqe_time,
                workflow.casci_time,
                workflow.casci_dmdm_time
            ]
        }
    ).to_csv(
        "compute_times.csv"
    )

    # plot the data
    x = np.linspace(0, 30, 1000)

    y = np.zeros_like(x)
    for i in range(args.num_states-1):
        y += gaussian(
            x,
            workflow._casci_results["exc_energies_ev"][i],
            workflow._casci_results["oscillator_strengths"][i]
        )

    plt.plot(x, y, label = "CASCI", alpha=0.6)

    y = np.zeros_like(x)
    for i in range(args.num_states-1):
        y += gaussian(
            x,
            workflow._casci_dmdm_results["exc_energies_ev"][i],
            workflow._casci_dmdm_results["oscillator_strengths"][i]
        )

    plt.plot(x, y, label = "CASCI + DMDM", alpha=0.6)

    y = np.zeros_like(x)
    for i in range(args.num_states-1):
        y += gaussian(
            x,
            workflow._vqe_results["exc_energies_ev"][i],
            workflow._vqe_results["oscillator_strengths"][i]
        )

    plt.plot(x, y, label = "VQE + DMDM (active space)", alpha=0.8, linestyle="--")

    y = np.zeros_like(x)
    for i in range(args.num_states-1):
        y += gaussian(
            x,
            workflow2._vqe_results["exc_energies_ev"][i],
            workflow2._vqe_results["oscillator_strengths"][i]
        )

    plt.plot(x, y, label = "VQE + DMDM", alpha=0.8, linestyle="dotted")



    plt.title(
        f"{args.molecule} uv-vis spectrum, {args.b}, CAS({args.num_active_orbitals},{args.num_active_electrons}), states: {args.num_states}"
    )

    plt.legend()
    plt.xlabel("Energy (eV)")
    plt.ylabel("Intensity (Oscillator Strength)")
    plt.savefig(f"{args.molecule}_uv-vis_spectrum_{args.b}_CAS({args.num_active_orbitals}_{args.num_active_electrons})_states_{args.num_states}.png")



if __name__ == "__main__":
    main()