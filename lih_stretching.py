"""
Stretching experiment of LiH
"""

import os
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
    gaussian
    )

from dmdm.interface import DMDM
import qrunch as qc

calculator = (qc.calculator_creator()
                .vqe()
                .iterative()
                .standard()
                .with_options(
                        options=qc.options.IterativeVqeOptions(max_iterations=500)
                )
                .choose_stopping_criterion(
                )
                .patience(
                        patience=10,
                        threshold=1e-10
                )
                .create()
)


vqe_dfs = []
casci_dfs = []
casci_dmdm_dfs = []

for factor in np.linspace(0.0, 4.0, 41):

    mol = [('Li', 0.0, 0.0, 0.0), ('H', 0.0, 0.0, 1.595+factor)]
    workflow = DMDMWorkflow(
        basis="cc-pvdz",
        molecule=mol,
        num_active_orbitals=2,
        num_active_electrons=2,
        num_states=4,
        mode=CalculationMode.BOTH,
        calculator=calculator,
        casci_like=True
    )

    # VQE only
    result_vqe = workflow.run_quantum_vqe()
    result_casci_dmdm = workflow.run_classical_casci_dmdm()
    result_casci = workflow.run_classical_casci()

    vqe_dfs.append(
        pd.DataFrame({
            "state": [*range(1, len(result_vqe["exc_energies_ev"])+1)],
            "energy_ev": result_vqe["exc_energies_ev"],
            "stretch": [factor] * len(result_vqe["exc_energies_ev"]),
            "oscillator_strength": result_vqe["oscillator_strengths"]
        })
    )
    
    casci_dfs.append(
        pd.DataFrame({
            "state": [*range(1, len(result_casci["exc_energies_ev"])+1)],
            "energy_ev": result_casci["exc_energies_ev"],
            "stretch": [factor] * len(result_casci["exc_energies_ev"]),
            "oscillator_strength": result_casci["oscillator_strengths"]
        })
    )

    casci__dmdm_dfs.append(
        pd.DataFrame({
            "state": [*range(1, len(result_casci_dmdm["exc_energies_ev"])+1)],
            "energy_ev": result_casci_dmdm["exc_energies_ev"],
            "stretch": [factor] * len(result_casci_dmdm["exc_energies_ev"]),
            "oscillator_strength": result_casci_dmdm["oscillator_strengths"]
        })
    )

vqe_df = pd.concat(vqe_dfs, axis=0).reset_index()
casci_df = pd.concat(casci_dfs, axis=0).reset_index()
casci_dmdm_df = pd.concat(casci_dmdm_dfs, axis=0).reset_index()
vqe_df.to_csv("lih_stretching/vqe_lih_stretching.csv")
casci_df.to_csv("lih_stretching/casci_lih_stretching.csv")
casci_dmdm_df.to_csv("lih_stretching/casci_dmdm_lih_stretching.csv")