""" Uitilities.
"""

import os
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
import copy
import re
import matplotlib.pyplot as plt
from typing import Any, Dict, List, Optional, Tuple
from enum import Enum
from pyscf import gto, scf, fci, mp, mcscf, ao2mo
import qrunch as qc
from qrunch.chemistry.reduced_density_matrices.reduced_density_matrix_calculator import ReducedDensityMatrixCalculator
from dmdm.interface import DMDM





def run_dalton(
    molecule_path : str,
    dalton_path : str,
    output_path : str,
):

    try:
        # Run a command (e.g., list files on Linux/Mac or dir on Windows)
        # Use a list for the command and arguments
        result = subprocess.run(
            [
                "dalton",
                "-dal" , dalton_path,
                "-mol", molecule_path,
                "-o", output_path
            ],
            capture_output=True,
            text=True,
            check=True
            )
        
        print("Return Code:", result.returncode)
        print("Output:\n", result.stdout)
        print("Errors:\n", result.stderr)
        
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}")
        print(f"Error output: {e.stderr}")
    except FileNotFoundError:
        print("Command not found!")



def parse_dalton_output(filename: str):
    """
    Parses Dalton output to extract excitation energies and oscillator strengths.
    Returns a pandas DataFrame.
    """
    try:
        with open(filename, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"Error: {filename} not found.")
        return None

    # Regex patterns
    # Matches the start of an excited state block
    state_pattern = r'@ Excited state no:\s+(\d+)\s+in symmetry\s+\d+\s+\(\s*\w+\s*\)\s+-\s+singlet excitation'
    
    # Matches excitation energy in atomic units (au)
    energy_pattern = r'@ Excitation energy\s*:\s+([\d\.E+-]+)\s+au'
    
    # Matches oscillator strength (LENGTH) for X, Y, Z
    # We capture the value and the polarization direction
    osc_pattern = r'@ Oscillator strength \(LENGTH\)\s*:\s+([\d\.E+-]+)\s+\((X|Y|Z)-polarization\)'

    # Find all state blocks
    states = []
    current_state = {}
    
    # Split content by state blocks to iterate
    # We use finditer to get the position of each state start
    state_matches = list(re.finditer(state_pattern, content))
    
    for i, match in enumerate(state_matches):
        state_num = int(match.group(1))
        
        # Define the start and end of this state's section
        start_idx = match.start()
        if i + 1 < len(state_matches):
            end_idx = state_matches[i+1].start()
        else:
            end_idx = len(content)
            
        section_text = content[start_idx:end_idx]
        
        # Extract Energy
        energy_match = re.search(energy_pattern, section_text)
        energy_au = float(energy_match.group(1)) if energy_match else None
        
        # Convert to eV (1 au = 27.2114 eV)
        energy_ev = energy_au * 27.2114 if energy_au else None
        
        # Extract Oscillator Strengths
        osc_matches = re.findall(osc_pattern, section_text)
        f_x = f_y = f_z = None
        
        for val, pol in osc_matches:
            val_float = float(val)
            if pol == 'X':
                f_x = val_float
            elif pol == 'Y':
                f_y = val_float
            elif pol == 'Z':
                f_z = val_float
        
        # Calculate Total Oscillator Strength (sum of components)
        f_total = (f_x or 0) + (f_y or 0) + (f_z or 0)

        states.append({
            'State': state_num,
            'Energy_au': energy_au,
            'Energy_eV': energy_ev,
            'Oscillator_F_X': f_x,
            'Oscillator_F_Y': f_y,
            'Oscillator_F_Z': f_z,
            'Oscillator_F_Total': f_total
        })

    if not states:
        print("No excited states found in the output file.")
        return None

    return pd.DataFrame(states)

def gaussian(x: np.array, energy: float, oscillator_strength: float, sigma: float=0.2):
    return oscillator_strength * np.exp(-(x - energy)**2 / (2 * sigma**2))


def scale_molecule(molecule: list, scale_factor: float, basis: str) -> gto.M:
    """
    Scales the atomic coordinates of a molecule by a given factor.
    
    Parameters:
    -----------
    molecule : list of tuples
        A list of tuples representing the atoms and their coordinates.
        Format: [('atom_name', x, y, z), ...]
    scale_factor : float
        The factor by which to stretch (greater than 1) or compress (less than 1) the molecule.
    basis: str
        The basis.

    Returns:
    --------
    gto.M
        A new molecule object with the scaled atomic coordinates.
    """
    # Extract the atom names and coordinates
    atom_names = [atom[0] for atom in molecule]
    coords = np.array([atom[1:] for atom in molecule])  # Coordinates (x, y, z)
    
    # Apply the scaling factor to all coordinates (stretching or compressing)
    new_coords = coords * scale_factor
    
    # Create a new molecule with the scaled coordinates
    new_molecule = gto.M(
        atom=[(atom_names[i], tuple(new_coords[i])) for i in range(len(molecule))],
        basis=basis  # Use a default basis set (you can change it if needed)
    )
    
    return new_molecule


def scale_molecule_v2(molecule: list, scale_factor: float) -> list:
    """
    Scales the atomic coordinates of a molecule by a given factor.
    
    Parameters:
    -----------
    molecule : list of lists or list of tuples
        A list of sequences representing the atoms and their coordinates.
        Format: [['atom_name', x, y, z], ...] or [('atom_name', x, y, z), ...]
    scale_factor : float
        The factor by which to stretch (greater than 1) or compress (less than 1) the molecule.

    Returns:
    --------
    list
        A new list of lists with the scaled atomic coordinates.
    """
    # We construct a new list of lists to ensure immutability of the input 
    # regardless of whether it was passed as tuples or lists.
    scaled_molecule = []
    
    for atom_data in molecule:
        # Extract the atom name (index 0)
        atom_name = atom_data[0]
        
        # Extract coordinates (indices 1 onwards) and scale them
        # Using a list comprehension for clarity and safety
        scaled_coords = [coord * scale_factor for coord in atom_data[1:]]
        
        # Append the new atom entry as a list
        scaled_molecule.append((atom_name,) + tuple(scaled_coords))
    
    return scaled_molecule


def write_dalton_molecule_file(molecule: list, filename: str, basis: str = "cc-pVDZ", 
                        molecule_name: str = None) -> None:
    """
    Writes molecule coordinates to a text file in the specified format.
    
    Parameters:
    -----------
    molecule : list of lists
        A list of lists representing the atoms and their coordinates.
        Format: [['atom_name', x, y, z], ...]
    filename : str
        Output file path
    basis : str
        Basis set name (default: cc-pVDZ)
    molecule_name : str, optional
        Name of the molecule (defaults to first atom element if not provided)
    """
    # Atomic charges (hardcoded for common elements - adjust as needed)
    atomic_charges = {
        'H': 1.0, 'He': 2.0, 'Li': 3.0, 'Be': 4.0, 'B': 5.0, 'C': 6.0,
        'N': 7.0, 'O': 8.0, 'F': 9.0, 'Ne': 10.0
    }
    
    # Group atoms by element type
    atom_groups = {}
    for atom in molecule:
        element = atom[0]
        if element not in atom_groups:
            atom_groups[element] = []
        atom_groups[element].append(atom[1:])  # Store only coordinates
    
    # Determine molecule name
    if molecule_name is None:
        # Build name from elements (simple approach)
        molecule_name = ''.join(atom_groups.keys())
    
    # Calculate total atoms
    total_atoms = sum(len(coords_list) for coords_list in atom_groups.values())
    
    # Write to file
    with open(filename, 'w') as f:
        # Header section
        f.write("BASIS\n")
        f.write(f"{basis}\n")
        f.write(f"{molecule_name}\n")
        f.write(f"using the {basis} basis\n")
        f.write(f"Atomtypes={len(atom_groups)} Nosymmetry\n")
        
        # Write each atom type
        for element, coords_list in atom_groups.items():
            charge = atomic_charges.get(element, 1.0)  # Default to 1.0 if unknown
            f.write(f"Charge={charge:.1f} Atoms={len(coords_list)}\n")
            
            for coords in coords_list:
                x, y, z = coords
                f.write(f"{element} {x:.10f} {y:.10f} {z:.10f}\n")


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
        scale_factor: float = 1.0,
        # VQE Specific Inputs
        vqe_max_iterations: int = 1000,
        vqe_patience: int = 10,
        vqe_threshold: float = 1e-10,
        verbose: int = 0,
        casci_like: bool = False,
        calculator: any = None,
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
            casci_like: Run DMDM like CASCI.
            calculator: qrunch calculator object for VQE
        """
        self.casci_like = casci_like
        self.basis = basis
        self.molecule_list = molecule
        self.num_active_orbitals = num_active_orbitals
        self.num_active_electrons = num_active_electrons 
        self.num_states = num_states
        self.mode = mode
        self.hartree_to_ev = 27.2114
        self.verbose = verbose
        self.calculator = calculator

        # VQE Config
        self.vqe_max_iterations = vqe_max_iterations
        self.vqe_patience = vqe_patience
        self.vqe_threshold = vqe_threshold

        # Results Storage
        self._casci_results: Optional[Dict] = None
        self._casscf_results: Optional[Dict] = None
        self._vqe_results: Optional[Dict] = None
        self._dmdm_casci: Optional[Any] = None
        self._dmdm_vqe: Optional[Any] = None
        
        # PEC Results Storage
        self._pec_casci: Optional[Dict] = None
        self._pec_casscf: Optional[Dict] = None
        self._pec_vqe: Optional[Dict] = None
        
        # Internal state flags
        self._casci_done = False
        self._casscf_done = False
        self._vqe_done = False
        self.scale_factor = scale_factor

        # Initial molecule setup
        self.molecule = scale_molecule(
            self.molecule_list,
            self.scale_factor,
            self.basis
        )


    def run_classical_casscf_average(self) -> Dict[str, Any]:
        """Run the classical CASSCF workflow using PySCF."""
        

        weights = np.ones(self.num_states)/self.num_states
        # 2. RHF & MP2
        mf = scf.RHF(self.molecule).run()
        mp2 = mp.MP2(mf).run()
        _, natorbs = mcscf.addons.make_natural_orbitals(mp2)

        # 3. CASSCF (for multiple roots, specify nroots > 1)
        mc = mcscf.CASSCF(
            mf,
            ncas=self.num_active_orbitals,
            nelecas=self.num_active_electrons
            ).state_average_(weights)
        mc.max_cycle = 100  # Increase the max cycles if needed
        mc.conv_tol = 1e-8  # Tighter convergence

        # Optionally, use natural orbitals as initial guess
        mc.mo_coeff = natorbs  # Uncomment if you want to use natural orbitals

        # 4. Run the CASSCF calculation to compute all roots
        e_tot, e_cas, ci, mo, mo_energy = mc.kernel()

        # 5. Integrals
        h_mo, _ = mc.get_h1eff()
        g_mo = mc.get_h2eff()
        g_mo = ao2mo.restore(1, g_mo, self.num_active_orbitals)

        # 6. RDM Reconstruction & DMDM
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
        E_core = e_tot - e_cas
        rdm_total_energies = rdm_active_energies + E_core

        # Dipole Integrals
        x_ao, y_ao, z_ao = self.molecule.intor('int1e_r', comp=3)
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

        # Store results
        self._casscf_results = {
            'exc_energies_ev': exc_energies[:self.num_states-1],
            'oscillator_strengths': osc_strengths[:self.num_states-1],
            'total_energies': rdm_total_energies,
            'e_cas': e_cas,
            'dmdm_obj': dmdm
        }
        self._casscf_done = True

        if self.verbose > 0:
            print("Done with CASSCF computations...\n")
        
        return self._casscf_results


    def run_classical_casscf(self) -> Dict[str, Any]:
        """Run the classical CASSCF workflow using PySCF."""

        weights = np.ones(self.num_states)/self.num_states
        # 2. RHF & MP2
        mf = scf.RHF(self.molecule).run()
        mp2 = mp.MP2(mf).run()
        _, natorbs = mcscf.addons.make_natural_orbitals(mp2)

        # 5. Integrals


        # 6. RDM Reconstruction & DMDM
        rdm_active_energies = []
        rdm_data_list = []
        energies_direct = []
        e_ground_state = 0.0

        for i in range(self.num_states):
            print("Optimising state: ", i)
            mc = mcscf.CASSCF(
            mf,
            ncas=self.num_active_orbitals,
            nelecas=self.num_active_electrons
            ).state_specific_(state=i)

            mc.max_cycle = 1000  # Increase the max cycles if needed
            mc.conv_tol = 1e-8  # Tighter convergence
            mc.mo_coeff = natorbs

            e_tot, e_cas, ci, mo, mo_energy = mc.kernel()
            if i == 0:
                e_ground_state = e_cas

            h_mo, _ = mc.get_h1eff()
            g_mo = mc.get_h2eff()
            g_mo = ao2mo.restore(1, g_mo, self.num_active_orbitals)

            ci_vec = ci
            rdm1, rdm2, rdm3, rdm4 = mc.fcisolver.make_rdm1234(ci_vec, self.num_active_orbitals, self.num_active_electrons)
            
            e1 = np.einsum('pq,pq', h_mo, rdm1)
            e2 = 0.5 * np.einsum('pqrs,pqrs', g_mo, rdm2)
            rdm_active_energies.append(e1 + e2)
            rdm_data_list.append((rdm1, rdm2, rdm3, rdm4))
            energies_direct.append(
                (e_cas - e_ground_state) * self.hartree_to_ev
            )

        rdm_active_energies = np.array(rdm_active_energies)
        E_core = e_tot - e_cas
        rdm_total_energies = rdm_active_energies + E_core

        # Dipole Integrals
        x_ao, y_ao, z_ao = self.molecule.intor('int1e_r', comp=3)
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

        # Store results
        self._casscf_results = {
            'exc_energies_direct_ev': np.array(energies_direct),
            'exc_energies_ev': exc_energies[:self.num_states-1],
            'oscillator_strengths': osc_strengths[:self.num_states-1],
            'total_energies': rdm_total_energies,
            'e_cas': e_cas,
            'dmdm_obj': dmdm
        }
        self._casscf_done = True

        if self.verbose > 0:
            print("Done with CASSCF computations...\n")
        
        return self._casscf_results

    # ==========================
    # CLASSICAL (CASCI) PATH
    # ==========================

    def run_classical_casci(self) -> Dict[str, Any]:
        """Run the classical CASCI workflow using PySCF."""

        # 2. RHF & MP2
        mf = scf.RHF(self.molecule).run()
        # mp2 = mp.MP2(mf).run()
        # _, natorbs = mcscf.addons.make_natural_orbitals(mp2)

        # 3. CASCI
        mc = mcscf.CASCI(mf, ncas=self.num_active_orbitals, nelecas=self.num_active_electrons)
        # mc = mcscf.CASSCF(mf, ncas=self.num_active_orbitals, nelecas=self.num_active_electrons)
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
        x_ao, y_ao, z_ao = self.molecule.intor('int1e_r', comp=3)
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
            'exc_energies_ev': exc_energies[:self.num_states-1],
            'oscillator_strengths': osc_strengths[:self.num_states-1],
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
        molecular_configuration = qc.build_molecular_configuration(molecule=self.molecule_list, basis_set=self.basis)
        
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
        # calculator = (
        #     qc.calculator_creator()
        #     .vqe()
        #     # .iterative()
        #     .iterative_with_orbital_optimization()
        #     .standard()
        #     .with_options(options=qc.options.IterativeVqeOptions(max_iterations=self.vqe_max_iterations))
        #     .choose_stopping_criterion()
        #     .patience(patience=self.vqe_patience, threshold=self.vqe_threshold)
        #     .create()
        # )

        # 6. Run Calculation
        result = self.calculator.calculate(ground_state_problem)

        # 7. Extract Integrals & RDMs
        self.molecule = molecular_configuration.pyscf_molecule
        mo_coeffs = problem_builder._restricted_builder._calculate_molecular_orbitals(molecular_configuration).alpha.coefficients

        # Transform AO to MO
        h_mo = one_electron_integral_transform(mo_coeffs, self.molecule.intor("int1e_kin") + self.molecule.intor("int1e_nuc"))
        g_mo = two_electron_integral_transform(mo_coeffs, self.molecule.intor("int2e"))

        # Overwrite active space with problem integrals
        cas_slice = slice(num_inactive_orbs, num_inactive_orbs + self.num_active_orbitals)
        h_mo[cas_slice, cas_slice] = ground_state_problem.electronic_structure_integrals.one_body_core_hamiltonian.alpha_alpha
        g_mo[cas_slice, cas_slice, cas_slice, cas_slice] = ground_state_problem.electronic_structure_integrals.two_body_electron_repulsion_integrals.alpha_alpha

        # 8. Calculate RDMs
        estimator = qc.estimator_creator().memory_restricted().with_precise_defaults().create()
        rdm_calculator = ReducedDensityMatrixCalculator(estimator=estimator)
        
        rdm1 = rdm_calculator.calculate_1_rdm(circuit=result.final_circuit, shots=None)
        rdm2 = rdm_calculator.calculate_2_rdm(circuit=result.final_circuit, shots=None)
        rdm3 = rdm_calculator.calculate_3_rdm(circuit=result.final_circuit, shots=None)
        rdm4 = rdm_calculator.calculate_4_rdm(circuit=result.final_circuit, shots=None)

        if self.casci_like == False:
            dmdm = DMDM(
                h_mo,
                g_mo,
                num_inactive_orbs,
                self.num_active_orbitals,
                num_virtual_orbs,
                molecular_configuration.number_of_electrons(),
                rdm1,
                rdm2=rdm2,
                rdm3=rdm3,
                rdm4=rdm4
            )
            coefs = mo_coeffs

        else:
            # 9. DMDM Initialization
            # Slice integrals for active space only (as per original script logic)
            h_mo_act = h_mo[cas_slice, cas_slice]
            g_mo_act = g_mo[cas_slice, cas_slice, cas_slice, cas_slice]

            dmdm = DMDM(
                h_mo_act,
                g_mo_act,
                0,
                self.num_active_orbitals, 
                0, 
                self.num_active_electrons,
                rdm1,
                rdm2=rdm2, 
                rdm3=rdm3,
                rdm4=rdm4
            )
            coefs = mo_coeffs[:, cas_slice]



        # 10. Dipole & Oscillator Strengths
        x, y, z = self.molecule.intor('int1e_r', comp=3)

        # Note: Original script used coefs = mo_coeffs[:, :num_active_electrons] which seems odd for spatial orbitals
        # Usually we slice by active orbitals. Using the active slice here for consistency.

        MO_DM = [
            one_electron_integral_transform(coefs, x),
            one_electron_integral_transform(coefs, y),
            one_electron_integral_transform(coefs, z)
        ]

        exc_energies_hartree = dmdm.get_excitation_energies()
        exc_energies_ev = exc_energies_hartree * self.hartree_to_ev
        osc_strengths = dmdm.get_oscillator_strength(MO_DM)

        self._vqe_results = {
            'exc_energies_ev': exc_energies_ev[:self.num_states-1],
            'oscillator_strengths': osc_strengths[:self.num_states-1],
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
        show_casscf: bool = True,
        show_vqe: bool = True,
        sigma: float = 0.2,
        title: Optional[str] = None,
        show: bool = True
    ) -> plt.Figure:
        """
        Plot the excitation spectrum.
        Can plot CASCI, VQE, or both for comparison.
        """
        if not show_casci and not show_casscf and not show_vqe:
            raise ValueError("Must show at least one method (CASCI, CASSCF or VQE).")

        x = np.linspace(0, 30, 1000)
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        legend_items = []

        # Plot CASCI
        if show_casci and self._casci_results:
            data = self._casci_results
            spectrum = np.zeros_like(x)
            for e, f in zip(data['exc_energies_ev'], data['oscillator_strengths']):
                spectrum += f * np.exp(-(x - e)**2 / (2 * sigma**2))
            
            ax.plot(x, spectrum, label="CASCI (Smooth)", color='green', alpha=0.8)
            ax.vlines(data['exc_energies_ev'], 0, data['oscillator_strengths'], color='blue', alpha=0.4, linewidth=1)
            legend_items.append("CASCI")
        
        if show_casscf and self._casscf_results:
            data = self._casscf_results
            spectrum = np.zeros_like(x)
            for e, f in zip(data['exc_energies_ev'], data['oscillator_strengths']):
                spectrum += f * np.exp(-(x - e)**2 / (2 * sigma**2))
            
            ax.plot(x, spectrum, label="CASSCF (Smooth)", color='blue', alpha=0.8)
            ax.vlines(data['exc_energies_ev'], 0, data['oscillator_strengths'], color='blue', alpha=0.4, linewidth=1)
            legend_items.append("CASSCF")

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

        if self.mode in [CalculationMode.CLASSICAL, CalculationMode.BOTH]:
            try:
                results['casscf'] = self.run_classical_casscf_average()
            except Exception as e:
                print(f"Error running CASSCF: {e}")

        if self.mode in [CalculationMode.QUANTUM, CalculationMode.BOTH]:
            try:
                results['vqe'] = self.run_quantum_vqe()
            except Exception as e:
                print(f"Error running VQE: {e}")

        if plot:
            self.plot_spectrum(show_casci='casci' in results, show_vqe='vqe' in results)

        return results

    def compute_pec(
        self,
        scale_range: Tuple[float, float],
        num_points: int = 10,
        method: str = "casci",
        plot: bool = False
        ) -> Dict[str, Any]:
        """
        Compute excited-state potential energy curves by scaling the molecule.
        
        Uses your scale_molecule function to generate displaced geometries.
        
        Args:
            scale_range: (min_scale, max_scale) scaling factors
            num_points: Number of geometry points to evaluate
            method: 'casci', 'casscf', 'casscf_average', or 'vqe'
            plot: Whether to plot the results immediately
        
        Returns:
            Dictionary containing scale_factors, energies for each state, and metadata
        """
        if method not in ['casci', 'casscf', 'casscf_average', 'vqe']:
            raise ValueError(f"Unknown method: {method}. Choose from 'casci', 'casscf', 'casscf_average', 'vqe'")
        
        if self.verbose > 0:
            print(f"Computing PEC for {method} method...")
        
        # Generate scale factors
        scale_factors = np.linspace(scale_range[0], scale_range[1], num_points)
        energies = np.zeros((num_points, self.num_states))
        distances = []  # Store actual distances if needed
        
        # Store original molecule
        original_molecule = self.molecule_list
        original_scale = self.scale_factor
        
        for i, scale in enumerate(scale_factors):
            if self.verbose > 1:
                print(f"  Computing point {i+1}/{num_points}: scale = {scale:.3f}")
            
            try:
                # Scale the molecule using your function
                self.molecule = scale_molecule(
                    original_molecule,
                    scale,
                    self.basis
                )
                
                # Update scale_factor for tracking
                self.scale_factor = scale
                
                # Run the chosen method
                if method == 'casci':
                    result = self.run_classical_casci()
                elif method == 'casscf':
                    result = self.run_classical_casscf()
                elif method == 'casscf_average':
                    result = self.run_classical_casscf_average()
                elif method == 'vqe':
                    result = self.run_quantum_vqe()
                
                # Extract excitation energies (first state is ground state at 0.0)
                exc_energies = result['exc_energies_ev']
                energies[i, :] = exc_energies
                
                # Store actual distance info if available (e.g., bond length)
                # This depends on your molecule - could calculate specific bond distances
                distances.append(scale)
                
            except Exception as e:
                if self.verbose > 0:
                    print(f"  Error at point {i} (scale={scale:.3f}): {e}")
                energies[i, :] = np.nan
        
        # Restore original molecule
        self.molecule = scale_molecule(original_molecule, original_scale, self.basis)
        self.scale_factor = original_scale
        
        # Store results
        pec_data = {
            'scale_factors': scale_factors,
            'energies': energies,
            'method': method,
            'scale_range': scale_range,
            'num_points': num_points,
            'states': list(range(self.num_states)),
            'molecule': original_molecule
        }
        
        # Store in appropriate attribute
        if method == 'casci':
            self._pec_casci = pec_data
        elif method in ['casscf', 'casscf_average']:
            self._pec_casscf = pec_data
        elif method == 'vqe':
            self._pec_vqe = pec_data
        
        if plot:
            self.plot_pec(method=method)
        
        return pec_data

    def plot_pec(
        self,
        method: Optional[str] = None,
        show_casci: bool = True,
        show_casscf: bool = True,
        show_vqe: bool = True,
        figsize: Tuple[int, int] = (10, 6),
        title: Optional[str] = None,
        show: bool = True
    ) -> plt.Figure:
        """
        Plot potential energy curves for one or multiple methods.
        
        Args:
            method: Specific method to plot ('casci', 'casscf', 'vqe')
            show_casci: Whether to show CASCI curves
            show_casscf: Whether to show CASSCF curves
            show_vqe: Whether to show VQE curves
            figsize: Figure size (width, height)
            title: Plot title
            show: Whether to display the plot
        
        Returns:
            matplotlib Figure object
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        colors = {
            'casci': 'green',
            'casscf': 'blue',
            'vqe': 'red'
        }
        
        linestyles = {
            'casci': '-',
            'casscf': '--',
            'vqe': '-.'
        }
        
        legend_labels = []
        
        # Helper to plot curves for a given method
        def plot_method_curves(pec_data, method_name, color, linestyle):
            if pec_data is None:
                return
            
            scale_factors = pec_data['scale_factors']
            energies = pec_data['energies']
            
            for state_idx in range(self.num_states):
                state_energies = energies[:, state_idx]
                
                # Handle NaN values
                valid_mask = ~np.isnan(state_energies)
                if np.any(valid_mask):
                    ax.plot(
                        scale_factors[valid_mask],
                        state_energies[valid_mask],
                        label=f"{method_name} State {state_idx}",
                        color=color,
                        linestyle=linestyle,
                        alpha=0.8,
                        linewidth=2
                    )
        
        # Plot each method if requested
        if method:
            # Plot specific method only
            if method == 'casci' and show_casci and self._pec_casci:
                plot_method_curves(self._pec_casci, "CASCI", colors['casci'], linestyles['casci'])
                legend_labels.append("CASCI")
            elif method == 'casscf' and show_casscf and self._pec_casscf:
                plot_method_curves(self._pec_casscf, "CASSCF", colors['casscf'], linestyles['casscf'])
                legend_labels.append("CASSCF")
            elif method == 'vqe' and show_vqe and self._pec_vqe:
                plot_method_curves(self._pec_vqe, "VQE", colors['vqe'], linestyles['vqe'])
                legend_labels.append("VQE")
        else:
            # Plot all requested methods
            if show_casci and self._pec_casci:
                plot_method_curves(self._pec_casci, "CASCI", colors['casci'], linestyles['casci'])
                legend_labels.append("CASCI")
            
            if show_casscf and self._pec_casscf:
                plot_method_curves(self._pec_casscf, "CASSCF", colors['casscf'], linestyles['casscf'])
                legend_labels.append("CASSCF")
            
            if show_vqe and self._pec_vqe:
                plot_method_curves(self._pec_vqe, "VQE", colors['vqe'], linestyles['vqe'])
                legend_labels.append("VQE")
        
        ax.set_xlabel("Scaling Factor", fontsize=12)
        ax.set_ylabel("Energy (eV)", fontsize=12)
        ax.set_title(
            title or f"Excited-State Potential Energy Curves ({', '.join(legend_labels)})",
            fontsize=14
        )
        ax.legend(loc='best', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        if show:
            plt.tight_layout()
            plt.show()
        
        return fig


    def compare_pec_methods(
        self,
        scale_range: Tuple[float, float],
        num_points: int = 10,
        methods: Optional[List[str]] = None,
        plot: bool = True
    ) -> Dict[str, Dict]:
        """
        Compute and compare PECs across multiple methods.
        
        Args:
            scale_range: (min_scale, max_scale) scaling factors
            num_points: Number of points
            methods: List of methods to compare (default: all available)
            plot: Whether to plot comparison
        
        Returns:
            Dictionary with PEC data for each method
        """
        if methods is None:
            methods = ['casci', 'casscf_average', 'vqe']
        
        results = {}
        
        for method in methods:
            try:
                if self.verbose > 0:
                    print(f"Computing PEC for {method}...")
                
                result = self.compute_pec(
                    scale_range=scale_range,
                    num_points=num_points,
                    method=method,
                    plot=False
                )
                results[method] = result
                
            except Exception as e:
                if self.verbose > 0:
                    print(f"Error computing PEC for {method}: {e}")
                results[method] = None
        
        if plot:
            self.plot_pec(show_casci='casci' in results, 
                         show_casscf=('casscf' in results) or ('casscf_average' in results), 
                         show_vqe='vqe' in results)
        
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