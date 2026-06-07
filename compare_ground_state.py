"""
Program to get VQE circuits and compare ground states
"""

import os
import argparse
from pathlib import Path
import pandas as pd
from pyscf import gto, scf, mcscf
from qchem.utils import  (
    MoleculeData,
    circuit_depth,
    get_all_gates,
)

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
        "-m",
        type=str,
        default="fast_vqe_as",
        help="VQE method name"
    )

    args = parser.parse_args()
    print(args)
    print("\n")

    # Create output folder
    Path(f"{args.o}").mkdir(parents=True, exist_ok=True)

    # Change the global working directory
    os.chdir(args.o)

    vqe_max_iterations = 500
    gate_selector = qc.gate_selector_creator().adapt().create()
    
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


    mol_name = args.molecule
    num_active_orbitals = args.num_active_orbitals
    num_active_electrons = args.num_active_electrons
    basis = args.b

    molecule = MoleculeData.molecules[mol_name]["coords"]

    mol = gto.M(
        atom=molecule,
        basis=basis
    )

    # -------------------------
    # 2. RHF
    # -------------------------
    mf = scf.RHF(mol).run()
    mc = mcscf.CASCI(mf, ncas=num_active_orbitals, nelecas=num_active_electrons)
    e_tot, _, _, _, _ = mc.kernel()

    print("CASCI e: ", e_tot)


    # VQE part
    molecular_configuration = qc.build_molecular_configuration(molecule=molecule, basis_set=basis)

    num_inactive_orbs = molecular_configuration.number_of_alpha_electrons() - (num_active_electrons // 2)
    num_spatial_orbs = molecular_configuration.number_of_spatial_orbitals()
    num_virtual_orbs = num_spatial_orbs - num_inactive_orbs - num_active_orbitals

    # 4. Build Problem
    problem_builder = (
        qc.problem_builder_creator()
        .ground_state()
        .standard()
        .add_problem_modifier()
        .active_space(
            number_of_active_spatial_orbitals=num_active_orbitals,
            number_of_active_alpha_electrons=num_active_electrons // 2
        )
        .create()
    )
    ground_state_problem = problem_builder.build_restricted(molecular_configuration)

    # 5. Setup Calculator (VQE)
    calculator = calculators[args.m]

    result = calculator.calculate(ground_state_problem)
    print("VQE e: ", result.total_energy.value)

    df = pd.DataFrame(
        {
            "molecule": [mol_name],
            "basis": [basis],
            "CAS": [f"CAS({num_active_electrons},{num_active_orbitals})"],
            "e_casci": [e_tot],
            "e_vqe": [result.total_energy.value],
            "vqe_method": [args.m],
            "final_circuit": [str(result.final_circuit)],
            "num_qubits": [result.final_circuit.num_qubits],
            # "num_gates_compact": [result.final_circuit.num_qubits],
            "num_gates_total": [len(get_all_gates(result.final_circuit))],
            "circuit_depth": [circuit_depth(result.final_circuit)]
        }
    ).to_csv(
        f"ground_state_{args.m}_{args.molecule}_{args.b}_CAS({args.num_active_electrons}_{args.num_active_orbitals}).csv"
    )




if __name__ == "__main__":
    main()