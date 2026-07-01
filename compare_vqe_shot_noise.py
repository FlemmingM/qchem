"""
Runs an experiment to compare VQE with different amounts of shot noise
"""

import argparse
import subprocess
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd
from qchem.utils import  (
    MoleculeData,
    DMDMWorkflow,
)

from dmdm.interface import DMDM
import qrunch as qc




def rmse(a: np.array, b: np.array) -> float:
    return np.sqrt(np.mean((a-b)**2))


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "molecule",
        type=str,
        help="molecule name"
    )

    parser.add_argument(
        "num_active_electrons",
        type=int,
        help="num_active_electrons"
    )

    parser.add_argument(
        "num_active_orbitals",
        type=int,
        help="num_active_orbitals"
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

    parser.add_argument(
        "-r",
        type=int,
        default=10,
        help="Number of repetitions"
    )

    args = parser.parse_args()

    print(args)

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
        .create()
    )

    workflow = DMDMWorkflow(
        basis=args.b,
        molecule=molecule,
        num_active_orbitals=args.num_active_orbitals,
        num_active_electrons=args.num_active_electrons,
        num_states=args.num_states,
        casci_like=True,
        calculator=calculator,
        shots=None
    )

    # take the time
    workflow.run_quantum_vqe(use_noisy_backend=False)

    pd.DataFrame(
        workflow._vqe_results
        ).to_csv(
            f"vqe_{args.molecule}_no_noise_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.csv"
        )

    shot_list = [1000, 5000, 10000, 50000, 100000]

    df_results = pd.DataFrame()
    for rep in range(args.r):
        for num_shots in shot_list:
            calculator = (
                qc.calculator_creator()
                .vqe()
                .iterative()
                .standard()
                .with_options(options=qc.options.IterativeVqeOptions(max_iterations=500))
                .create()
            )

            workflow = DMDMWorkflow(
                basis=args.b,
                molecule=molecule,
                num_active_orbitals=args.num_active_orbitals,
                num_active_electrons=args.num_active_electrons,
                num_states=args.num_states,
                casci_like=True,
                calculator=calculator,
                shots=num_shots
            )

            # take the time
            workflow.run_quantum_vqe(
                use_noisy_backend=True,
                )

            # compare to no shots
            # Add reps and shots to results
            workflow._vqe_results["shots"] = [num_shots] * len(workflow._vqe_results["exc_energies_ev"])
            workflow._vqe_results["rep"] = [rep] * len(workflow._vqe_results["exc_energies_ev"])

            # append dataframe
            df_results = pd.concat(
                [
                    df_results,
                    pd.DataFrame(workflow._vqe_results),
                ],
                axis=0
            )

    df_results.reset_index().drop("index", axis=1).to_csv(
        f"vqe_{args.molecule}_{args.b}_shot_noise_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.csv"
    )

if __name__ == "__main__":
    main()