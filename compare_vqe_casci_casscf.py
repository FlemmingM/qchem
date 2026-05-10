"""
Runs an experiment to compare VQE with CACI, CASSCF
"""

import argparse
import subprocess
import os
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
    run_dalton,
    )

from dmdm.interface import DMDM
import qrunch as qc


# qc.register_license_file("/home/flemming/Nextcloud/Cherimoya/training/master_cs/ms_project/code/qchem/license_fm.txt")




# Path("my_folder").mkdir(parents=True, exist_ok=True)


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
        "-dal",
        type=str,
        help="Dalton file path"
    )

    parser.add_argument(
        "-mol",
        type=str,
        help="Molecule file path for dalton program"
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
    dalton_out = args.mol.split("/")[-1].replace(".mol", "") + ".out"

    # Create output folder
    Path(f"{args.o}").mkdir(parents=True, exist_ok=True)

    # Change the global working directory
    os.chdir(args.o)

    # Run CASSCF
    run_dalton(
        "../"+args.mol,
        "../"+args.dal,
        dalton_out
    )


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
        num_active_orbitals=2,
        num_active_electrons=2,
        num_states=4,
        mode=CalculationMode.BOTH,
        vqe_max_iterations=50000,
        casci_like=True,
        calculator=calculator
    )

    workflow.run_quantum_vqe()
    workflow.run_classical_casci()


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
        num_active_orbitals=2,
        num_active_electrons=2,
        num_states=4,
        mode=CalculationMode.BOTH,
        vqe_max_iterations=50000,
        casci_like=False,
        calculator=calculator
    )

    workflow2.run_quantum_vqe()

    # process CASSCF
    df = parse_dalton_output(dalton_out)


    # plot the data
    x = np.linspace(0, 30, 1000)
    

    y = np.zeros_like(x)
    for i in range(len(df)):
        y += gaussian(
            x,
            df.loc[i, "Energy_eV"],
            df.loc[i, "Oscillator_F_Total"]
        )
    plt.plot(x, y, label = "CASSCF", alpha=0.6)

    y = np.zeros_like(x)
    for i in range(len(workflow._casci_results["exc_energies_ev"])):
        y += gaussian(
            x,
            workflow._casci_results["exc_energies_ev"][i],
            workflow._casci_results["oscillator_strengths"][i]
        )

    plt.plot(x, y, label = "CASCI", alpha=0.6)

    y = np.zeros_like(x)
    for i in range(len(workflow._vqe_results["exc_energies_ev"])):
        y += gaussian(
            x,
            workflow._vqe_results["exc_energies_ev"][i],
            workflow._vqe_results["oscillator_strengths"][i]
        )

    plt.plot(x, y, label = "VQE + DMDM (active space)", alpha=0.8, linestyle="--")

    y = np.zeros_like(x)
    for i in range(len(workflow2._vqe_results["exc_energies_ev"])):
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
    # plt.show()

if __name__ == "__main__":
    main()