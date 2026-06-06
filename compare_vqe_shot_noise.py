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
from qiskit_aer.noise import NoiseModel, depolarizing_error




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

    print("Running VQE vs. CASCI comparison\n")
    print(args)
    print("\n")

    # Create output folder
    Path(f"{args.o}").mkdir(parents=True, exist_ok=True)

    # Change the global working directory
    os.chdir(args.o)


    # Create a noise model
    noise_model = NoiseModel()
    noise_model.add_all_qubit_quantum_error(depolarizing_error(0.01, 1), ['sx'])
    noise_model.add_all_qubit_quantum_error(depolarizing_error(0.02, 2), ['cx'])


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
        casci_like=True,
        calculator=calculator,
        shots=None
    )

    # take the time
    workflow.run_quantum_vqe(use_noisy_backend=False)
    print(f"VQE (cas-like) time (seconds): {workflow.vqe_time}")

    rdms_true = workflow.vqe_rdms

    pd.DataFrame(
        workflow._vqe_results
        ).to_csv(
            f"vqe_{args.molecule}_no_noise_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.csv"
        )

    # shot_list = [10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000]

    shot_list = [1000, 5000]

    rdm1_errors = []
    rdm2_errors = []

    df_results = pd.DataFrame()
    for rep in range(args.r):
        for num_shots in shot_list:
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
                casci_like=True,
                calculator=calculator,
                shots=num_shots
            )

            # take the time
            workflow.run_quantum_vqe(
                use_noisy_backend=True,
                noise_model=noise_model
                )

            # compare to no shots
            current_rdms = workflow.vqe_rdms
            rdm1_error = rmse(rdms_true[0], current_rdms[0])
            rdm2_error = rmse(rdms_true[1], current_rdms[1])

            print(f"Shots: {num_shots}\n")
            print(f"RDM1 error: {rdm1_error}")
            print(f"RDM2 error: {rdm2_error}")

            rdm1_errors.append(rdm1_error)
            rdm2_errors.append(rdm2_error)

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
        
    # save the data
    pd.DataFrame(
        {
            "shots": shot_list * args.r,
            "rdm1_errors": rdm1_errors,
            "rdm2_errors": rdm2_errors
        }
    ).to_csv(
        f"{args.molecule}_rdm12_errors_{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.csv"
    )

    df_results.reset_index().drop("index", axis=1).to_csv(
        f"vqe_{args.molecule}_{args.b}_shot_noise_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.csv"
    )




if __name__ == "__main__":
    main()