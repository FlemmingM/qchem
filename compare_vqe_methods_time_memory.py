"""
Runs an experiment to get spectral properties, time and memory consumption.
"""

import argparse
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
        "-method",
        type=str,
        default="casci",
        help="Method"
    )

    parser.add_argument(
        "-tm",
        # type=bool,
        action='store_true',
        default=False,
        help="Save times and memory t/f"
    )

    parser.add_argument(
        "-s",
        # type=bool,
        action='store_false',
        default=True,
        help="Save spectra t/f"
    )

    args = parser.parse_args()

    # Create output folder
    Path(f"{args.o}").mkdir(parents=True, exist_ok=True)

    # Change the global working directory
    os.chdir(args.o)


    vqe_max_iterations = 500


    # Get molecule
    molecule = MoleculeData.molecules[args.molecule]["coords"]
    method = args.method

    # Define the calculators


    if method in ["fast_vqe_as", "casci", "casci_dmdm"]:
        calculator = (qc.calculator_creator()
                .vqe()
                .iterative()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .create()
        )

    if method == "oofast_vqe_as":
        calculator = (qc.calculator_creator()
                .vqe()
                .iterative_with_orbital_optimization()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .create()
        )

    if method == "adapt_vqe_as":
        # Gate selector
        gate_selector = qc.gate_selector_creator().adapt().create()
        calculator = (qc.calculator_creator()
                .vqe()
                .iterative()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .with_gate_selector(gate_selector)
                .create()
        )

    if method == "ooadapt_vqe_as": 
        # Gate selector
        gate_selector = qc.gate_selector_creator().adapt().create()
        calculator = (qc.calculator_creator()
                .vqe()
                .iterative_with_orbital_optimization()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .with_gate_selector(gate_selector)
                .create()
        )

    if method == "fast_vqe": 
        calculator = (qc.calculator_creator()
                .vqe()
                .iterative()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .create()
        )

    if method == "oofast_vqe": 
        calculator = (
                qc.calculator_creator()
                .vqe()
                .iterative_with_orbital_optimization()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .create()
        )

    if method == "adapt_vqe": 
        # Gate selector
        gate_selector = qc.gate_selector_creator().adapt().create()
        calculator = (
                qc.calculator_creator()
                .vqe()
                .iterative()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .with_gate_selector(gate_selector)
                .create()
        )

    if method == "ooadapt_vqe":
        # Gate selector
        gate_selector = qc.gate_selector_creator().adapt().create()
        calculator = (
                qc.calculator_creator()
                .vqe()
                .iterative_with_orbital_optimization()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .with_gate_selector(gate_selector)
                .create()
        )

    methods = []
    times = []
    times_method = []
    memory_footprints_total = []
    memory_footprints_method = []
    memory_footprints_dmdm = []


    if "_as" in method:
        casci_like = True
    else:
        casci_like = False

    workflow = DMDMWorkflow(
        basis=args.b,
        molecule=molecule,
        num_active_orbitals=args.num_active_orbitals,
        num_active_electrons=args.num_active_electrons,
        num_states=args.num_states,
        casci_like=casci_like,
        calculator=calculator,
    )

    if method == "casci":
        print("running CASCI")
        workflow.run_classical_casci()
        memory_footprints_total.append(workflow.mem_total)
        memory_footprints_method.append(workflow.mem_method)
        memory_footprints_dmdm.append(None)
        methods.append("casci")
        times.append(workflow.casci_time)
        times_method.append(workflow.casci_time)

        if args.s:
            pd.DataFrame(
                    workflow._casci_results
                ).dropna().to_csv(
                    f"casci_{args.molecule}_{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.csv"
                )

    if method == "casci_dmdm":
        workflow.run_classical_casci_dmdm()
        memory_footprints_total.append(workflow.mem_total)
        memory_footprints_method.append(workflow.mem_method)
        memory_footprints_dmdm.append(workflow.mem_dmdm)
        methods.append("casci_dmdm")
        times.append(workflow.casci_dmdm_time)
        times_method.append(workflow.casci_dmdm_time_method)
        
        if args.s:
            pd.DataFrame(
                    workflow._casci_dmdm_results
                ).dropna().to_csv(
                    f"casci_dmdm_{args.molecule}_{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.csv"
                )

    if "vqe" in method:
        # Run VQE method
        workflow.run_quantum_vqe()
        memory_footprints_total.append(workflow.mem_total)
        memory_footprints_method.append(workflow.mem_method)
        memory_footprints_dmdm.append(workflow.mem_dmdm)


        # Add times
        methods.append(f"{method}")
        times.append(workflow.vqe_time)
        times_method.append(workflow.vqe_time_method)
        if args.s:
            pd.DataFrame(
                    workflow._vqe_results
                ).dropna().to_csv(
                    f"{method}_{args.molecule}_{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.csv"
                )

    if args.tm:
        pd.DataFrame(
            {
                "method": methods,
                "time": times,
                "time_method": times_method,
                "memory_dmdm_mb" : memory_footprints_dmdm,
                "memory_method_mb": memory_footprints_method,
                "memory_total_mb": memory_footprints_total,
                "basis": [f"{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}"] * len(methods),
                "molecule": [args.molecule] * len(methods)
            }
        ).to_csv(
            f"compute_times_memory_{method}_{args.molecule}_{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.csv"
        )


if __name__ == "__main__":
    main()