""" Uitilities.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Any, Dict, List, Optional, Tuple
from enum import Enum
from pyscf import gto, scf, fci, mp, mcscf, ao2mo
import qrunch as qc
from qrunch.chemistry.reduced_density_matrices.reduced_density_matrix_calculator import ReducedDensityMatrixCalculator
from dmdm.interface import DMDM

class MoleculeData:
    molecules = {
    "H2O": {
        "coords": [
            ("O", 0.0, 0.0, 0.1035174918),
            ("H", 0.0, 0.7955612117, -0.4640237459),
            ("H", 0.0, -0.7955612117, -0.4640237459),
        ],
        "minimal": {
            "num_active_orbs": 4,
            "num_active_alpha_electrons": 2
        },
        "larger": {
            "num_active_orbs": 6,
            "num_active_alpha_electrons": 3
        }
    },
    "LiH": {
        "coords": [
            ("Li", 0.0, 0.0, 0.0),
            ("H", 0.0, 0.0, 1.595),
        ],
        "minimal": {
            "num_active_orbs": 2,
            "num_active_alpha_electrons": 1
        },
        "larger": {
            "num_active_orbs": 4,
            "num_active_alpha_electrons": 2
        },
        "full": {
            "num_active_orbs": 6,
            "num_active_alpha_electrons": 2
        }
    },
    "HF": {
        "coords": [
            ("H", 0.0, 0.0, 0.0),
            ("F", 0.0, 0.0, 0.917),
        ],
        "minimal": {
            "num_active_orbs": 2,
            "num_active_alpha_electrons": 1
        },
        "larger": {
            "num_active_orbs": 6,
            "num_active_alpha_electrons": 3
        }
    },
    "N2": {
        "coords": [
            ("N", 0.0, 0.0, -0.550),
            ("N", 0.0, 0.0, 0.550),
        ],
        "minimal": {
            "num_active_orbs": 6,
            "num_active_alpha_electrons": 3
        },
        "larger": {
            "num_active_orbs": 10,
            "num_active_alpha_electrons": 5
        }
    },
    "CO": {
        "coords": [
            ("C", 0.0, 0.0, 0.0),
            ("O", 0.0, 0.0, 1.128),
        ],
        "minimal": {
            "num_active_orbs": 4,
            "num_active_alpha_electrons": 2
        },
        "larger": {
            "num_active_orbs": 8,
            "num_active_alpha_electrons": 4
        }
    },
    "CH4": {
        "coords": [
            ("C", 0.0, 0.0, 0.0),
            ("H", 0.0, 0.0, 1.085),
            ("H", 0.0, 0.943, -0.362),
            ("H", 0.817, -0.471, -0.362),
            ("H", -0.817, -0.471, -0.362),
        ],
        "minimal": {
            "num_active_orbs": 4,
            "num_active_alpha_electrons": 2
        },
        "larger": {
            "num_active_orbs": 8,
            "num_active_alpha_electrons": 4
        }
    },
    "NH3": {
        "coords": [
            ("N", 0.0, 0.0, 0.114),
            ("H", 0.0, 0.938, -0.342),
            ("H", 0.812, -0.469, -0.342),
            ("H", -0.812, -0.469, -0.342),
        ],
        "minimal": {
            "num_active_orbs": 4,
            "num_active_alpha_electrons": 2
        },
        "larger": {
            "num_active_orbs": 8,
            "num_active_alpha_electrons": 4
        }
    },
    "BeH2": {
        "coords": [
            ("Be", 0.0, 0.0, 0.0),
            ("H", 0.0, 0.0, 1.330),
            ("H", 0.0, 0.0, -1.330),
        ],
        "minimal": {
            "num_active_orbs": 2,
            "num_active_alpha_electrons": 1
        },
        "larger": {
            "num_active_orbs": 4,
            "num_active_alpha_electrons": 2
        }
    },
    "F2": {
        "coords": [
            ("F", 0.0, 0.0, -0.700),
            ("F", 0.0, 0.0, 0.700),
        ],
        "minimal": {
            "num_active_orbs": 4,
            "num_active_alpha_electrons": 2
        },
        "larger": {
            "num_active_orbs": 8,
            "num_active_alpha_electrons": 4
        }
    },
}


class CalculationMode(Enum):
    CLASSICAL = "classical"
    QUANTUM = "quantum"
    BOTH = "both"


class DMDMWorkflow:
    """
    Unified workflow for computing excitation spectra using either:
    1. Classical CASCI (PySCF)
    2. Quantum VQE (qrunch/qchem)
    
    Both paths feed into the same DMDM analysis and plotting routines.
    """

    def __init__(
        self,
        basis: str,
        molecule: list[tuple],
        num_active_orbitals: int,
        num_active_electrons: int,
        num_states: int,
        mode: CalculationMode = CalculationMode.BOTH,
        # VQE Specific Inputs
        vqe_max_iterations: int = 1000,
        vqe_patience: int = 10,
        vqe_threshold: float = 1e-10,
        verbose: int = 0
    ):
        """
        Initialize the workflow.
        
        Args:
            basis: Basis set string (e.g., 'aug-cc-pVDZ').
            molecule: molecule coordinates as list of tuples.
            num_active_orbitals: Number of active orbitals.
            num_active_electrons: Number of active electrons.
            num_states: Number of roots/states.
            mode: Whether to run Classical, Quantum, or Both.
            vqe_max_iterations: Max iterations for VQE optimizer.
            vqe_patience: Patience for VQE stopping criterion.
            vqe_threshold: Threshold for VQE convergence.
            verbose: Verbosity level.
        """
        self.basis = basis
        self.molecule = molecule
        self.num_active_orbitals = num_active_orbitals
        self.num_active_electrons = num_active_electrons
        self.num_states = num_states
        self.mode = mode
        self.hartree_to_ev = 27.2114
        self.verbose = verbose

        # VQE Config
        self.vqe_max_iterations = vqe_max_iterations
        self.vqe_patience = vqe_patience
        self.vqe_threshold = vqe_threshold

        # Results Storage
        self._casci_results: Optional[Dict] = None
        self._vqe_results: Optional[Dict] = None
        self._dmdm_casci: Optional[DMDM] = None
        self._dmdm_vqe: Optional[DMDM] = None
        
        # Internal state flags
        self._casci_done = False
        self._vqe_done = False

    # ==========================
    # CLASSICAL (CASCI) PATH
    # ==========================

    def run_classical_casci(self) -> Dict[str, Any]:
        """Run the classical CASCI workflow using PySCF."""
                

        mol = gto.M(
            atom=self.molecule,
            basis=self.basis
        )

        # 2. RHF & MP2
        mf = scf.RHF(mol).run()
        mp2 = mp.MP2(mf).run()
        # _, natorbs = mcscf.addons.make_natural_orbitals(mp2)

        # 3. CASCI
        mc = mcscf.CASCI(mf, ncas=self.num_active_orbitals, nelecas=self.num_active_electrons)
        mc.fcisolver.nroots = self.num_states
        # mc.mo_coeff = natorbs # Optional: use natural orbitals
        
        e_tot, e_cas, ci, mo, mo_energy = mc.kernel()
        
        # 4. Integrals
        h_mo, _ = mc.get_h1eff()
        g_mo = mc.get_h2eff()
        g_mo = ao2mo.restore(1, g_mo, self.num_active_orbitals)

        # 5. RDM Reconstruction & DMDM
        rdm_active_energies = []
        rdm_data_list = []

        for i in range(self.num_states):
            ci_vec = ci[i]
            rdm1, rdm2, rdm3, rdm4 = mc.fcisolver.make_rdm1234(ci_vec, self.num_active_orbitals, self.num_active_electrons)
            
            e1 = np.einsum('pq,pq', h_mo, rdm1)
            e2 = 0.5 * np.einsum('pqrs,pqrs', g_mo, rdm2)
            rdm_active_energies.append(e1 + e2)
            rdm_data_list.append((rdm1, rdm2, rdm3, rdm4))

        rdm_active_energies = np.array(rdm_active_energies)
        E_core = e_tot[0] - e_cas[0]
        rdm_total_energies = rdm_active_energies + E_core


        # Dipole Integrals
        x_ao, y_ao, z_ao = mol.intor('int1e_r', comp=3)
        cas_slice = slice(mc.ncore, mc.ncore + self.num_active_orbitals)
        C_cas = mo[:, cas_slice]
        
        x_cas = one_electron_integral_transform(C_cas, x_ao)
        y_cas = one_electron_integral_transform(C_cas, y_ao)
        z_cas = one_electron_integral_transform(C_cas, z_ao)
        MO_DM = [x_cas, y_cas, z_cas]

        # Initialize DMDM
        gs_rdm = rdm_data_list[0]
        dmdm = DMDM(
            h_mo,
            g_mo,
            0,
            self.num_active_orbitals,
            0,
            self.num_active_electrons,
            gs_rdm[0],
            rdm2=gs_rdm[1],
            rdm3=gs_rdm[2],
            rdm4=gs_rdm[3]
        )

        # Get energies
        exc_energies = dmdm.get_excitation_energies() * self.hartree_to_ev

        # Get Oscillator Strengths
        osc_strengths = dmdm.get_oscillator_strength(MO_DM)

        self._casci_results = {
            'exc_energies_ev': np.concatenate(([0.0], exc_energies[:self.num_states-1])),
            'oscillator_strengths': np.concatenate(([0.0], osc_strengths[:self.num_states-1])),
            'total_energies': rdm_total_energies,
            'e_cas': e_cas,
            'dmdm_obj': dmdm
        }
        self._casci_done = True
        if self.verbose > 0:
            print("Done with CASCI computations...\n")
        return self._casci_results

    # ==========================
    # QUANTUM (VQE) PATH
    # ==========================

    def run_quantum_vqe(self) -> Dict[str, Any]:
        """Run the VQE workflow using qrunch/qchem."""

        # 3. Build Configuration
        molecular_configuration = qc.build_molecular_configuration(molecule=self.molecule, basis_set=self.basis)
        
        num_inactive_orbs = molecular_configuration.number_of_alpha_electrons() - (self.num_active_electrons // 2)
        num_spatial_orbs = molecular_configuration.number_of_spatial_orbitals()
        num_virtual_orbs = num_spatial_orbs - num_inactive_orbs - self.num_active_orbitals

        print(f"  Inactive: {num_inactive_orbs}, Active: {self.num_active_orbitals}, Virtual: {num_virtual_orbs}")

        # 4. Build Problem
        problem_builder = (
            qc.problem_builder_creator()
            .ground_state()
            .standard()
            .add_problem_modifier()
            .active_space(
                number_of_active_spatial_orbitals=self.num_active_orbitals,
                number_of_active_alpha_electrons=self.num_active_electrons // 2
            )
            .create()
        )
        ground_state_problem = problem_builder.build_restricted(molecular_configuration)

        # 5. Setup Calculator (VQE)
        gate_selector = qc.gate_selector_creator().adapt().create()
        calculator = (
            qc.calculator_creator()
            .vqe()
            .iterative()
            .standard()
            .with_options(options=qc.options.IterativeVqeOptions(max_iterations=self.vqe_max_iterations))
            .choose_stopping_criterion()
            .patience(patience=self.vqe_patience, threshold=self.vqe_threshold)
            .create()
        )

        # 6. Run Calculation
        result = calculator.calculate(ground_state_problem)

        # 7. Extract Integrals & RDMs
        mol = molecular_configuration.pyscf_molecule
        mo_coeffs = problem_builder._restricted_builder._calculate_molecular_orbitals(molecular_configuration).alpha.coefficients

        # Transform AO to MO
        h_mo = one_electron_integral_transform(mo_coeffs, mol.intor("int1e_kin") + mol.intor("int1e_nuc"))
        g_mo = two_electron_integral_transform(mo_coeffs, mol.intor("int2e"))

        # Overwrite active space with problem integrals
        cas_slice = slice(num_inactive_orbs, num_inactive_orbs + self.num_active_orbitals)
        h_mo[cas_slice, cas_slice] = ground_state_problem.electronic_structure_integrals.one_body_core_hamiltonian.alpha_alpha
        # g_mo[cas_slice, cas_slice, cas_slice, cas_slice] = ground_state_problem.electronic_structure_integrals.two_body_electron_repulsion_integrals.alpha_alpha # Simplified slicing logic
        # Note: The original script did 4D slicing. Ensure shapes match.
        g_mo[cas_slice, cas_slice, cas_slice, cas_slice] = ground_state_problem.electronic_structure_integrals.two_body_electron_repulsion_integrals.alpha_alpha

        # 8. Calculate RDMs
        estimator = qc.estimator_creator().memory_restricted().with_precise_defaults().create()
        rdm_calculator = ReducedDensityMatrixCalculator(estimator=estimator)
        
        rdm1 = rdm_calculator.calculate_1_rdm(circuit=result.final_circuit, shots=None)
        rdm2 = rdm_calculator.calculate_2_rdm(circuit=result.final_circuit, shots=None)
        rdm3 = rdm_calculator.calculate_3_rdm(circuit=result.final_circuit, shots=None)
        rdm4 = rdm_calculator.calculate_4_rdm(circuit=result.final_circuit, shots=None)

        # 9. DMDM Initialization
        # Slice integrals for active space only (as per original script logic)
        h_mo_act = h_mo[cas_slice, cas_slice]
        g_mo_act = g_mo[cas_slice, cas_slice, cas_slice, cas_slice]

        dmdm = DMDM(
            h_mo_act, g_mo_act, 0, self.num_active_orbitals, 0, self.num_active_electrons,
            rdm1, rdm2=rdm2, rdm3=rdm3, rdm4=rdm4
        )

        # 10. Dipole & Oscillator Strengths
        x, y, z = mol.intor('int1e_r', comp=3)
        # Note: Original script used coefs = mo_coeffs[:, :num_active_electrons] which seems odd for spatial orbitals
        # Usually we slice by active orbitals. Using the active slice here for consistency.
        coefs = mo_coeffs[:, cas_slice] 
        
        MO_DM = [
            one_electron_integral_transform(coefs, x),
            one_electron_integral_transform(coefs, y),
            one_electron_integral_transform(coefs, z)
        ]

        exc_energies_hartree = dmdm.get_excitation_energies()
        exc_energies_ev = exc_energies_hartree * self.hartree_to_ev
        osc_strengths = dmdm.get_oscillator_strength(MO_DM)

        self._vqe_results = {
            'exc_energies_ev': np.concatenate(([0.0], exc_energies_ev[:self.num_states-1])),
            'oscillator_strengths': np.concatenate(([0.0], osc_strengths[:self.num_states-1])),
            'dmdm_obj': dmdm
        }
        self._vqe_done = True
        if self.verbose > 0:
            print("Done with VQE computations...\n")
        return self._vqe_results

    # ==========================
    # PLOTTING & ANALYSIS
    # ==========================

    def plot_spectrum(
        self,
        show_casci: bool = True,
        show_vqe: bool = True,
        sigma: float = 0.2,
        title: Optional[str] = None,
        show: bool = True
    ) -> plt.Figure:
        """
        Plot the excitation spectrum.
        Can plot CASCI, VQE, or both for comparison.
        """
        if not show_casci and not show_vqe:
            raise ValueError("Must show at least one method (CASCI or VQE).")

        x = np.linspace(0, 30, 1000)
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        legend_items = []

        # Plot CASCI
        if show_casci and self._casci_results:
            data = self._casci_results
            spectrum = np.zeros_like(x)
            for e, f in zip(data['exc_energies_ev'], data['oscillator_strengths']):
                spectrum += f * np.exp(-(x - e)**2 / (2 * sigma**2))
            
            ax.plot(x, spectrum, label="CASCI (Smooth)", color='blue', alpha=0.8)
            ax.vlines(data['exc_energies_ev'], 0, data['oscillator_strengths'], color='blue', alpha=0.4, linewidth=1)
            legend_items.append("CASCI")

        # Plot VQE
        if show_vqe and self._vqe_results:
            data = self._vqe_results
            spectrum = np.zeros_like(x)
            for e, f in zip(data['exc_energies_ev'], data['oscillator_strengths']):
                spectrum += f * np.exp(-(x - e)**2 / (2 * sigma**2))
            
            ax.plot(x, spectrum, label="VQE (Smooth)", color='red', linestyle='--', alpha=0.8)
            ax.vlines(data['exc_energies_ev'], 0, data['oscillator_strengths'], color='red', alpha=0.4, linewidth=1, linestyle='--')
            legend_items.append("VQE")

        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel("Intensity (Oscillator Strength)")
        ax.set_title(title or f"Excitation Spectrum: {', '.join(legend_items)}")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 30)

        if show:
            plt.show()
        
        return fig

    def run_comparison(self, plot: bool = True) -> Dict[str, Any]:
        """
        Convenience method to run both workflows (if enabled) and plot the comparison.
        """
        results = {}

        if self.mode in [CalculationMode.CLASSICAL, CalculationMode.BOTH]:
            try:
                results['casci'] = self.run_classical_casci()
            except Exception as e:
                print(f"Error running CASCI: {e}")

        if self.mode in [CalculationMode.QUANTUM, CalculationMode.BOTH]:
            try:
                results['vqe'] = self.run_quantum_vqe()
            except Exception as e:
                print(f"Error running VQE: {e}")

        if plot:
            self.plot_spectrum(show_casci='casci' in results, show_vqe='vqe' in results)

        return results


def get_hf_gse_from_mol(
    molecule: list,
    basis: str
    ) -> float:
    """
    Function to compute ground the state energy using FCI.
    """
    mol = gto.M(
    atom=molecule,
    basis=basis,
    unit="Angstrom",
    )

    hf_energy = mol.RHF().run()
    moller_plesset = mp.MP2(hf_energy).run()
    noons, natorbs = mcscf.addons.make_natural_orbitals(moller_plesset)
    cisolver = fci.FCI(mol, natorbs)
    fci_gse, fcivec = cisolver.kernel()

    return fci_gse, fcivec


def one_electron_integral_transform(c: np.ndarray, int1e: np.ndarray) -> np.ndarray:
    """
    Transform one-electron integrals from the atomic orbital (AO) basis to a
    molecular orbital (MO) basis using a coefficient matrix.

    Parameters
    ----------
    c : np.ndarray
        Coefficient matrix of shape ``(n_ao, n_mo)`` that expands molecular
        orbitals in terms of atomic orbitals. Each column corresponds to a
        molecular orbital expressed as a linear combination of AO basis
        functions.

    int1e : np.ndarray
        One-electron integral tensor in the AO basis with shape ``(n_ao, n_ao)``.
        Typical examples include the kinetic-energy matrix or nuclear-attraction
        matrix.

    Returns
    -------
    np.ndarray
        The transformed one-electron integral matrix in the MO basis,
        with shape ``(n_mo, n_mo)``. The operation performed is::

            int1e_MO[i, j] = Σ_a Σ_b C[a, i] * C[b, j] * int1e[a, b]

    Notes
    -----
    The function uses ``np.einsum`` with an explicit contraction path for
    efficiency. The path ``[(0, 2), (0, 1)]`` first contracts the first
    coefficient matrix with the integral tensor, then contracts the second
    coefficient matrix.

    Examples
    --------
    >>> C = np.random.rand(7, 5)          # 7 AOs → 5 MOs
    >>> int1e = np.random.rand(7, 7)     # AO one-electron integrals
    >>> int1e_mo = one_electron_integral_transform(C, int1e)
    >>> int1e_mo.shape
    (5, 5)
    """
    return np.einsum(
        "ai,bj,ab->ij",
        c,
        c,
        int1e,
        optimize=["einsum_path", (0, 2), (0, 1)],
    )


def two_electron_integral_transform(c: np.ndarray, int2e: np.ndarray) -> np.ndarray:
    """
    Transform two-electron integrals from the atomic orbital (AO) basis to a
    molecular orbital (MO) basis using a coefficient matrix.

    Parameters
    ----------
    c : np.ndarray
        Coefficient matrix of shape ``(n_ao, n_mo)`` that expands molecular
        orbitals in terms of atomic orbitals. Each column corresponds to a
        molecular orbital expressed as a linear combination of AO basis
        functions.

    int2e : np.ndarray
        Two-electron integral tensor in the AO basis with shape
        ``(n_ao, n_ao, n_ao, n_ao)``. This tensor encodes electron-repulsion
        integrals ⟨ab|cd⟩ in the AO basis.

    Returns
    -------
    np.ndarray
        The transformed two-electron integral tensor in the MO basis,
        with shape ``(n_mo, n_mo, n_mo, n_mo)``. The transformation follows::

            int2e_MO[i, j, k, l] =
                Σ_a Σ_b Σ_c Σ_d C[a, i] * C[b, j] *
                C[c, k] * C[d, l] * int2e[a, b, c, d]

    Notes
    -----
    ``np.einsum`` is used with a manually-specified contraction order for
    optimal performance. The path ``[(0, 4), (0, 3), (0, 2), (0, 1)]`` contracts
    each coefficient matrix sequentially with the four-index integral tensor.

    Examples
    --------
    >>> C = np.random.rand(7, 5)           # 7 AOs → 5 MOs
    >>> int2e = np.random.rand(7, 7, 7, 7) # AO two-electron integrals
    >>> int2e_mo = two_electron_integral_transform(C, int2e)
    >>> int2e_mo.shape
    (5, 5, 5, 5)
    """
    return np.einsum(
        "ai,bj,ck,dl,abcd->ijkl",
        c,
        c,
        c,
        c,
        int2e,
        optimize=[
            "einsum_path",
            (0, 4),
            (0, 3),
            (0, 2),
            (0, 1),
        ],
    )