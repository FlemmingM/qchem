import pyscf
import numpy as np
import qrunch as qc

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
    # Define the molecule and basis set.
    molecule = [
        ("O", 0.0, 0.0, 0.1035174918),
        ("H", 0.0, 0.7955612117, -0.4640237459),
        ("H", 0.0, -0.7955612117, -0.4640237459),
    ]

    basis = "sto-3g"

    # Build the molecular configuration
    molecular_configuration = qc.build_molecular_configuration(
        molecule=molecule,
        basis_set=basis,
    )

    # Set active space parameters
    # (4, 4) active space:
    # num_active_orbs = 4
    # num_active_alpha_electrons = 2
    num_active_orbs = 4
    num_active_alpha_electrons = 2

    num_inactive_orbs = molecular_configuration.number_of_alpha_electrons() - num_active_alpha_electrons
    num_virtual_orbs = molecular_configuration.number_of_spatial_orbitals() - num_inactive_orbs - num_active_orbs
    num_spatial_orbs = molecular_configuration.number_of_spatial_orbitals()

    print("Number of inactive orbitals: ", num_inactive_orbs)
    print("Number of active orbitals:   ", num_active_orbs)
    print("Number of virtual orbitals:  ", num_virtual_orbs)
    print("Number of spatial orbitals:  ", num_spatial_orbs)

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
    fast_vqe_calculator = (
        qc.calculator_creator()
        .vqe()
        .iterative()
        .standard()
        .create()
    )

    oofast_vqe_calculator = (
        qc.calculator_creator()
        .vqe()
        .iterative_with_orbital_optimization()
        .standard()
        .create()
    )

    gate_selector = (
        qc.gate_selector_creator()
        .adapt()
        .create()
    )

    adapt_vqe_calculator = (
        qc.calculator_creator()
        .vqe()
        .iterative()
        .standard()
        .with_gate_selector(gate_selector)
        .create()
    )

    ooadapt_vqe_calculator = (
        qc.calculator_creator()
        .vqe()
        .iterative_with_orbital_optimization()
        .standard()
        .with_gate_selector(gate_selector)
        .create()
    )

    calculator = oofast_vqe_calculator
    result = calculator.calculate(ground_state_problem)

    # Get the PySCF molecule and MO coefficients
    # Note that the MO coefficients are for the full space, so we will need to transform the integrals
    # and then overwrite the active space block with the integrals from the ground_state_problem.
    mol = molecular_configuration.pyscf_molecule
    mo_coeffs = problem_builder._restricted_builder._calculate_molecular_orbitals(molecular_configuration).alpha.coefficients

    # Transform the one- and two-electron AO integrals to the MO basis.
    h_mo = one_electron_integral_transform(mo_coeffs, mol.intor("int1e_kin") + mol.intor("int1e_nuc"))
    g_mo = two_electron_integral_transform(mo_coeffs, mol.intor("int2e"))


    # Get the MO integrals from the problem builder in the AS.
    h_mo[
        num_inactive_orbs:num_inactive_orbs+num_active_orbs,
        num_inactive_orbs:num_inactive_orbs+num_active_orbs
    ] = (
        ground_state_problem.electronic_structure_integrals.one_body_core_hamiltonian.alpha_alpha
        # + one_electron_integral_transform(mo_coeffs, mol.intor("int1e_nuc"))[
        #     num_inactive_orbs:num_inactive_orbs+num_active_orbs,
        #     num_inactive_orbs:num_inactive_orbs+num_active_orbs
        # ]
    )
    g_mo[
        num_inactive_orbs:num_inactive_orbs+num_active_orbs,
        num_inactive_orbs:num_inactive_orbs+num_active_orbs,
        num_inactive_orbs:num_inactive_orbs+num_active_orbs,
        num_inactive_orbs:num_inactive_orbs+num_active_orbs
    ] = ground_state_problem.electronic_structure_integrals.two_body_electron_repulsion_integrals.alpha_alpha

    print("One-body integrals shape:    ", h_mo.shape)
    print("Two-body integrals shape:    ", g_mo.shape)

    # Create the Estimator, requires a PRO license.
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
        molecular_configuration.number_of_electrons(),
        rdm1,
        rdm2=rdm2,
        rdm3=rdm3,
        rdm4=rdm4,
    )

    # Calculate excitation energies and oscillator strengths.
    x, y, z = mol.intor('int1e_r', comp=3)
    MO_DM = [one_electron_integral_transform(mo_coeffs, x), one_electron_integral_transform(mo_coeffs, y), one_electron_integral_transform(mo_coeffs, z)]

    exc = dmdm.get_excitation_energies()
    osc = dmdm.get_oscillator_strength(MO_DM)

    print(exc)
    print(osc)

if __name__ == "__main__":
    main()
