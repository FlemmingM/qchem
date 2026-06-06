"""
Runs an experiment to compare VQE with CACI, CASSCF
"""

import argparse
import subprocess
import os
import tracemalloc
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
    build_spectrum,
    spectral_similarity,
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
        "-mp2",
        type=bool,
        default=False,
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


    vqe_max_iterations = 500


    # Get molecule
    molecule = MoleculeData.molecules[args.molecule]["coords"]


    # Gate selector
    gate_selector = qc.gate_selector_creator().adapt().create()


    # Define the calculators

    calculators = {
        "fast_vqe_as":
                qc.calculator_creator()
                .vqe()
                .iterative()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .choose_stopping_criterion(
                )
                .patience(
                        patience=10,
                        threshold=1e-10
                )
                # .choose_minimizer()
                # .scipy(method="L-BFGS-B", options=qc.options.ScipyMinimizerOptions(jacobian_step_size=1e-9))
                .create(),

        "oofast_vqe_as":
                qc.calculator_creator()
                .vqe()
                .iterative_with_orbital_optimization()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .create(),

        "adapt_vqe_as":
                qc.calculator_creator()
                .vqe()
                .iterative()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .with_gate_selector(gate_selector)
                .create(),

        "ooadapt_vqe_as": 
                qc.calculator_creator()
                .vqe()
                .iterative_with_orbital_optimization()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .with_gate_selector(gate_selector)
                .create(),

        "fast_vqe":
                qc.calculator_creator()
                .vqe()
                .iterative()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .choose_stopping_criterion(
                )
                .patience(
                        patience=10,
                        threshold=1e-10
                )
                # .choose_minimizer()
                # .scipy(method="L-BFGS-B", options=qc.options.ScipyMinimizerOptions(jacobian_step_size=1e-9))
                .create(),

        "oofast_vqe":
                qc.calculator_creator()
                .vqe()
                .iterative_with_orbital_optimization()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .create(),

        "adapt_vqe":
                qc.calculator_creator()
                .vqe()
                .iterative()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .with_gate_selector(gate_selector)
                .create(),

        "ooadapt_vqe": 
                qc.calculator_creator()
                .vqe()
                .iterative_with_orbital_optimization()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=vqe_max_iterations,)
                )
                .with_gate_selector(gate_selector)
                .create()
    }


    methods = []
    times = []
    times_method = []
    memory_footprints_total = []
    memory_footprints_method = []

    # for the plot
    x = np.linspace(0, 40, 1000)
    spectral_similarities = []

    for name, calculator in calculators.items():
        if "_as" in name:
            casci_like=True
        else:
            casci_like=False

        workflow = DMDMWorkflow(
            basis=args.b,
            molecule=molecule,
            num_active_orbitals=args.num_active_orbitals,
            num_active_electrons=args.num_active_electrons,
            num_states=args.num_states,
            mode=CalculationMode.BOTH,
            casci_like=casci_like,
            calculator=calculator,
            mp2=args.mp2
        )

        if name == "fast_vqe_as":

            # run classical reference
            workflow.run_classical_casci()
            memory_footprints_total.append(workflow.mem_total)
            memory_footprints_method.append(workflow.mem_method)

            # add stats and save file
            pd.DataFrame(
                workflow._casci_results
            ).to_csv(
                f"casci_{args.molecule}_{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.csv"
            )

            methods.append("casci")
            times.append(workflow.casci_time)
            times_method.append(workflow.casci_time)

            # Plot
            spectrum_casci = build_spectrum(
                workflow._casci_results["exc_energies_ev"],
                workflow._casci_results["oscillator_strengths"],
                x
                )
            spectral_similarities.append(spectral_similarity(
                spectrum_casci, spectrum_casci, x
            ))
            plt.plot(x, spectrum_casci, label = "CASCI", alpha=0.6, color="black", linestyle="--")


            workflow.run_classical_casci_dmdm()
            memory_footprints_total.append(workflow.mem_total)
            memory_footprints_method.append(workflow.mem_method)

            # add stats and save file
            pd.DataFrame(
                workflow._casci_dmdm_results
            ).to_csv(
                f"casci_dmdm_{args.molecule}_{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.csv"
            )

            methods.append("casci_dmdm")
            times.append(workflow.casci_dmdm_time)
            times_method.append(workflow.casci_dmdm_time_method)

            # Plot
            spectrum_casci_dmdm = build_spectrum(
                workflow._casci_dmdm_results["exc_energies_ev"],
                workflow._casci_dmdm_results["oscillator_strengths"],
                x
                )
            spectral_similarities.append(spectral_similarity(
                spectrum_casci, spectrum_casci_dmdm, x
            ))
            plt.plot(x, spectrum_casci_dmdm, label = "CASCI + DMDM", alpha=0.6)

        # Run VQE method
        workflow.run_quantum_vqe()
        memory_footprints_total.append(workflow.mem_total)
        memory_footprints_method.append(workflow.mem_method)

        # Save dataframe

        try:
            pd.DataFrame(
                workflow._vqe_results
            ).to_csv(
                f"{name}_{args.molecule}_{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.csv"
            )
        except:
            import json
            with open(f"{name}_{args.molecule}_{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.json", 'w') as f:
                json.dump(data, f)

        # Add times
        methods.append(f"{name}")
        times.append(workflow.vqe_time)
        times_method.append(workflow.vqe_time_method)

        with open(f"compute_times_memory_{name}_{args.molecule}_{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.txt", "w") as outfile:
            outfile.write(f"method: {name}\n")
            outfile.write(f"time: {workflow.vqe_time}\n")
            outfile.write(f"time_method: {workflow.vqe_time_method}\n")
            outfile.write(f"memory_total_mb: {workflow.mem_total}\n")
            outfile.write(f"memory_method_mb: {workflow.mem_method}\n")

        plot_name = name.replace("_", " ").upper()

        # Plot
        spectrum = build_spectrum(
                workflow._vqe_results["exc_energies_ev"],
                workflow._vqe_results["oscillator_strengths"],
                x
            )
        spectral_similarities.append(spectral_similarity(
                spectrum, spectrum_casci, x
            ))

        plt.plot(x, spectrum, label = f"{plot_name}", alpha=0.6)


    # Finalise the plot
    plt.title(
        f"VQE Comparison{args.molecule} uv-vis spectrum, {args.b}, CAS({args.num_active_orbitals},{args.num_active_electrons}), states: {args.num_states}"
    )

    plt.legend()
    plt.xlabel("Energy (eV)")
    plt.ylabel("Intensity (Oscillator Strength)")
    plt.savefig(f"{args.molecule}_uv-vis_spectrum_{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.png")


    # save the compute times and memory footprints
    pd.DataFrame(
        {
            "method": methods,
            "time": times,
            "time_method": times_method,
            "memory_total_mb" : memory_footprints_total,
            "memory_method_mb": memory_footprints_method,
            "spectral_similarity": spectral_similarities,
            "basis": [f"{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}"] * len(methods),
            "molecule": [args.molecule] * len(methods)
        }
    ).to_csv(
        f"compute_times_memory_{args.molecule}_{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals})_states_{args.num_states}.csv"
    )


if __name__ == "__main__":
    main()