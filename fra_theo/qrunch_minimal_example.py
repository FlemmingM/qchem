# import os
# import pyscf
import numpy as np
import qrunch as qc
# from pyscf import mcscf, mp, fci

from qrunch.chemistry.reduced_density_matrices.reduced_density_matrix_calculator import (
    ReducedDensityMatrixCalculator,
)
from dmdm.interface import DMDM

def one_electron_integral_transform(C: np.ndarray, int1e: np.ndarray) -> np.ndarray:
    return np.einsum("ai,bj,ab->ij", C, C, int1e, optimize=["einsum_path", (0, 2), (0, 1)])

def two_electron_integral_transform(C: np.ndarray, int2e: np.ndarray) -> np.ndarray:
    return np.einsum(
        "ai,bj,ck,dl,abcd->ijkl", C, C, C, C, int2e, optimize=["einsum_path", (0, 4), (0, 3), (0, 2), (0, 1)]
    )

def main() -> None:

    qc.setup_logger(qc.LOGGER_ESSENTIAL)
    print(f"Qrunch version: {qc.__version__}")

    molecule = [
        ("Li", 0.0, 0.0, 0.0),
        ("H", 0.0, 0.0, 1.67),
    ]
    basis = "cc-pVDZ"

    # Build the molecular configuration
    molecular_configuration = qc.build_molecular_configuration(
        molecule=molecule,
        basis_set=basis,
        units="angstrom"
    )

    # Set active space parameters
    # (2, 2) active space:
    num_active_orbs = 2
    num_active_alpha_electrons = 1
    # Full space
    # num_active_orbs = molecular_configuration.number_of_spatial_orbitals()
    # num_active_alpha_electrons = molecular_configuration.number_of_alpha_electrons()

    num_inactive_orbs = molecular_configuration.number_of_alpha_electrons() - num_active_alpha_electrons
    num_virtual_orbs = molecular_configuration.number_of_spatial_orbitals() - num_inactive_orbs - num_active_orbs
    num_spatial_orbs = molecular_configuration.number_of_spatial_orbitals()
    num_electrons = molecular_configuration.number_of_electrons()

    print(f"Number of inactive orbitals:  {num_inactive_orbs}")
    print(f"Number of active orbitals:    {num_active_orbs}")
    print(f"Number of virtual orbitals:   {num_virtual_orbs}")
    print(f"Number of spatial orbitals:   {num_spatial_orbs}")
    print(f"Active space:                ({num_active_orbs}, {2*num_active_alpha_electrons})")

    # Build the ground state problem.
    problem_builder = (
        qc.problem_builder_creator()
        .ground_state()
        .standard()
        .add_problem_modifier()
        .active_space(
            number_of_active_spatial_orbitals=num_active_orbs,
            number_of_active_alpha_electrons=num_active_alpha_electrons,
        )
        .create()
    )
    ground_state_problem = problem_builder.build_restricted(molecular_configuration)

    # To test different calculators.

    fci_calculator = (
        qc.calculator_creator()
        .configuration_interaction()
        .standard()
        .create()
    )

    # os.makedirs(f"LiH/{basis}/{2*num_active_alpha_electrons}_{num_active_orbs}", exist_ok=True)
    fast_vqe_calculator = (
        qc.calculator_creator()
        .vqe()
        .iterative_with_orbital_optimization()
        # .iterative()
        .standard()
        .with_options(
            options=qc.options.IterativeVqeOptions(max_iterations=1000,)
        )
        # .choose_minimizer()
        # .scipy(method="L-BFGS-B", options=qc.options.ScipyMinimizerOptions(jacobian_step_size=1e-9, function_tolerance=1e-12, gradient_tolerance=1e-12, max_function_evaluations=1000,))
        .choose_stopping_criterion()
        .patience(
            patience=10,
            threshold=1e-10,
        )
        # .choose_data_persister_manager()
        # .file_persister(f"LiH/{basis}/{2*num_active_alpha_electrons}_{num_active_orbs}")
        .create()
    )

    calculator = fast_vqe_calculator
    result = calculator.calculate(ground_state_problem)
    print(f"E(Qrunch) = {result.total_energy}")
    fci_results = fci_calculator.calculate(ground_state_problem)
    print(f"E(FCI [Qrunch]) = {fci_results.total_energy}")
    print(f"Energy difference: {result.total_energy - fci_results.total_energy}")

    # Get the PySCF molecule, MO coefficients and integrals
    mol = molecular_configuration.pyscf_molecule
    mo_coeffs = problem_builder._restricted_builder._calculate_molecular_orbitals(molecular_configuration).alpha.coefficients

    h_mo = one_electron_integral_transform(mo_coeffs, mol.intor("int1e_kin") + mol.intor("int1e_nuc"))
    g_mo = two_electron_integral_transform(mo_coeffs, mol.intor("int2e"))

    print("MO coefficients shape:       ", mo_coeffs.shape)
    print("One-body integrals shape:    ", h_mo.shape)
    print("Two-body integrals shape:    ", g_mo.shape)

    # Create the Estimator that require a PRO license.
    estimator = (
        qc.estimator_creator()
        .memory_restricted()
        .with_precise_defaults()
        .create()
    )

    # Create the RDM calculator and calculate the RDMs.
    rdm_calculator = ReducedDensityMatrixCalculator(
        estimator=estimator,
    )
    rdm1 = rdm_calculator.calculate_1_rdm(circuit=result.final_circuit, shots=None)
    rdm2 = rdm_calculator.calculate_2_rdm(circuit=result.final_circuit, shots=None)
    rdm3 = rdm_calculator.calculate_3_rdm(circuit=result.final_circuit, shots=None)
    rdm4 = rdm_calculator.calculate_4_rdm(circuit=result.final_circuit, shots=None)

    # Initialize DMDM
    dmdm = DMDM(
        h_mo,
        g_mo,
        num_inactive_orbs,
        num_active_orbs,
        num_virtual_orbs,
        num_electrons,
        rdm1,
        rdm2=rdm2,
        rdm3=rdm3,
        rdm4=rdm4,
    )

    # Calculate excitation energies and oscillator strengths.
    x, y, z = mol.intor('int1e_r', comp=3)
    MO_DM = [one_electron_integral_transform(mo_coeffs, x), one_electron_integral_transform(mo_coeffs, y), one_electron_integral_transform(mo_coeffs, z)]

    exc = dmdm.get_excitation_energies()
    uvvis_output = dmdm.get_uvvis_output(MO_DM)
    print(uvvis_output)

if __name__ == "__main__":
    main()
