""" Uitilities.
"""

import os
import time
from pathlib import Path
import psutil
import subprocess
import numpy as np
import pandas as pd
import copy
import re
import matplotlib.pyplot as plt
from typing import Any, Dict, List, Optional, Tuple
from enum import Enum
from pyscf import gto, scf, fci, mp, mcscf, ao2mo
import pyscf
import qrunch as qc
from qrunch.chemistry.reduced_density_matrices.reduced_density_matrix_calculator import ReducedDensityMatrixCalculator
from dmdm.interface import DMDM
from typing import Dict, Any, Optional
from math import comb
from collections import defaultdict


def add_missed_transitions(df_casci, df_vqe, penalty=10.0):
    """
    Add penalty-shifted entries for transitions present in CASCI but missing from VQE.
    
    For each transition in df_casci whose index is not referenced in df_vqe['ref_idx'],
    a new row is appended to df_vqe with pred_E = ref_E + penalty.
    
    Parameters:
    - df_casci: Reference DataFrame (CASCI), indexed by transition
    - df_vqe: Predicted DataFrame (VQE), must contain 'ref_idx', 'ref_E', 'pred_E'
    - penalty: Energy shift (eV) added to ref_E for missing transitions
    
    Returns:
    - df_vqe with missing transitions appended
    """
    matched_idxs = set(df_vqe["ref_idx"].tolist())
    
    missed_rows = []
    for idx in df_casci.index:
        if idx not in matched_idxs:
            ref_e = df_casci.loc[idx, "exc_energies_ev"]
            ref_f = df_casci.loc[idx, "oscillator_strengths"]
            missed_rows.append({
                "ref_idx": idx,
                "ref_E": ref_e,
                "pred_E": ref_e + penalty,
                "ref_f": ref_f,
                "pred_f": 0.0
            })
    
    if missed_rows:
        df_vqe = pd.concat([df_vqe, pd.DataFrame(missed_rows)], axis=0, ignore_index=True)
    
    return df_vqe



# helper function to build spectra and comparison tables
def compare_spectra(df_true, dfs_pred, labels, sigma=0.2, tol=2.0, k=10, x_range=(0, 80), x_points=1000, col: str="oscillator_strengths"):
    """
    Compare predicted spectra against a reference using Weighted RMSE and Spectral Similarity.
    
    Args:
        df_true (pd.DataFrame): Reference dataframe with 'exc_energies_ev' and 'oscillator_strengths'.
        dfs_pred (list of pd.DataFrame): List of predicted dataframes.
        labels (list of str): Method labels (first entry should be the reference method).
        sigma (float): Gaussian broadening width in eV for spectral similarity.
        tol (float): Maximum energy difference (eV) for greedy matching in RMSE.
        k (int): Number of top-K reference states for RMSE.
        x_range (tuple): Energy range (min, max) for spectral construction.
        x_points (int): Number of points in the energy grid.
        col (str): The column used as y axis.
        
    Returns:
        tuple:
            - results_df (pd.DataFrame): Columns 'method', 'rmse', 'spectral_similarity'.
            - x_grid (np.ndarray): The energy grid used for spectra.
            - y_ref (np.ndarray): The reference spectrum array.
            - y_preds (dict): Dictionary mapping method labels to their spectrum arrays.
    """
    x = np.linspace(x_range[0], x_range[1], x_points)
    
    # Build reference spectrum once
    y_ref = build_spectrum(
        df_true["exc_energies_ev"],
        df_true[col],
        x,
        sigma=sigma
    )
    
    # Initialize lists for results
    rmses = [0.0]
    similarities = [1.0]
    y_preds = {labels[0]: y_ref}  # Store reference spectrum under its label
    
    for i, df in enumerate(dfs_pred):
        method_label = labels[i + 1] # Skip the first label which is for the reference
        
        
        # 2. Spectral Similarity & Store Spectrum
        y_pred = build_spectrum(
            df["exc_energies_ev"],
            df[col],
            x,
            sigma=sigma
        )
        
        sim = spectral_similarity(y_ref, y_pred, x)
        similarities.append(sim)
        
        # Store the spectrum
        y_preds[method_label] = y_pred
    
    results_df = pd.DataFrame({
        "method": labels,
        # "rmse": rmses,
        "spectral_similarity": similarities
    })
    
    return results_df, x, y_ref, y_preds


def match_spectral_peaks(
    df_ref: pd.DataFrame,
    dfs_pred: list,
    top_k: int = 10,
    tol_eV: float = 0.2,
    penalty: float = 10.0
) -> tuple[list[pd.DataFrame], list[pd.DataFrame]]:
    """
    Matches the top k strongest peaks from a reference spectrum against multiple
    predicted spectra using a greedy energy-proximity algorithm.

    Args:
        df_ref: DataFrame with columns ['exc_energies_ev', 'oscillator_strengths', 'rotational_strengths'].
        dfs_pred: List of DataFrames (one per method) with same columns.
        top_k: Number of strongest reference peaks to attempt to match.
        tol_eV: Maximum energy difference (eV) to consider a peak a match.

    Returns:
        Tuple of two lists:
        - matched_dfs: One DataFrame per method, one row per reference peak.
        - unmatched_dfs: One DataFrame per method, containing predicted peaks
          that were NOT matched to any reference peak.
    """

    # --- 1. Select top-K reference peaks by oscillator strength ---
    df_ref_sorted = (
        df_ref.sort_values("oscillator_strengths", ascending=False)
        .head(top_k)
        .reset_index(drop=True)
        .dropna(subset=["exc_energies_ev"])
    )

    n_targets = len(df_ref_sorted)
    if n_targets == 0:
        return [], []

    matched_dfs = []
    unmatched_dfs = []

    for idx, df_pred in enumerate(dfs_pred):
        df_pred_sorted = df_pred.sort_values("exc_energies_ev").reset_index(drop=True)
        used_pred_indices = set()

        # --- 2. Build rows greedily ---
        rows = []
        for i, ref_row in df_ref_sorted.iterrows():
            ref_e = ref_row["exc_energies_ev"]
            ref_f = ref_row["oscillator_strengths"]
            ref_r = ref_row["rotational_strengths"]

            best_dist = np.inf
            best_j = -1

            for j, pred_row in df_pred_sorted.iterrows():
                if j in used_pred_indices:
                    continue
                dist = abs(pred_row["exc_energies_ev"] - ref_e)
                if dist < best_dist:
                    best_dist = dist
                    best_j = j

            if best_dist <= tol_eV and best_j != -1:
                # --- Match found ---
                used_pred_indices.add(best_j)
                pred_row = df_pred_sorted.iloc[best_j]
                rows.append({
                    "ref_exc_energies_ev": ref_e,
                    "ref_oscillator_strengths": ref_f,
                    "ref_rotational_strengths": ref_r,
                    "pred_exc_energies_ev": pred_row["exc_energies_ev"],
                    "pred_oscillator_strengths": pred_row["oscillator_strengths"],
                    "pred_rotational_strengths": pred_row["rotational_strengths"],
                    "delta_e": best_dist,
                    "delta_o": abs(ref_f - pred_row["oscillator_strengths"]),
                    "matched": True,
                })
            else:
                # --- No match → NaN ---
                rows.append({
                    "ref_exc_energies_ev": ref_e,
                    "ref_oscillator_strengths": ref_f,
                    "ref_rotational_strengths": ref_r,
                    "pred_exc_energies_ev": np.nan,
                    "pred_oscillator_strengths": np.nan,
                    "pred_rotational_strengths": np.nan,
                    "delta_e": penalty,
                    "delta_o": 0.0,
                    "matched": False,
                })

        result_df = pd.DataFrame(rows).sort_values("ref_exc_energies_ev").reset_index(drop=True)

        # --- 3. Collect unmatched predicted peaks ---
        unmatched_mask = ~df_pred_sorted.index.isin(used_pred_indices)
        unmatched_df = (
            df_pred_sorted.loc[unmatched_mask, ["exc_energies_ev", "oscillator_strengths", "rotational_strengths"]]
            .sort_values("exc_energies_ev")
            .reset_index(drop=True)
        )

        # --- 4. Attach summary stats ---
        n_matched = result_df["matched"].sum()
        result_df.attrs = {
            "method_idx": idx,
            "n_targets": n_targets,
            "n_matched": n_matched,
            "match_ratio": n_matched / n_targets,
            "avg_energy_error_eV": result_df["delta_e"].mean(),
            "avg_osc_error": (
                result_df.loc[result_df["matched"], "ref_oscillator_strengths"]
                - result_df.loc[result_df["matched"], "pred_oscillator_strengths"]
            ).abs().mean(),
        }
        matched_dfs.append(result_df)
        unmatched_dfs.append(unmatched_df)

    return matched_dfs, unmatched_dfs


def match_spectral_peaks_v2(
    df_ref: pd.DataFrame, 
    dfs_pred: list, 
    top_k: int = 10, 
    tol_eV: float = 0.2, 
    unmatched_penalty: float = 10.0,
) -> tuple:
    """
    Matches the top k strongest peaks from a reference spectrum against multiple 
    predicted spectra using a greedy energy-proximity algorithm.
    
    Args:
        df_ref: DataFrame with columns ['exc_energies_ev', 'oscillator_strengths'].
        dfs_pred: List of DataFrames (one per method) with same columns.
        top_k: Number of strongest reference peaks to attempt to match.
        tol_eV: Maximum energy difference (eV) to consider a peak a match.
        unmatched_penalty: Value used for error calculation if a peak cannot be matched.
        
    Returns:
        A tuple containing:
        - results_list: List of dicts, one per predicted DF, containing matched errors and stats.
        - global_matches: Details of all matches (useful for debugging/plotting).
    """
    
    # 1. Identify the Top K Reference Peaks (sorted by oscillator strength)
    df_ref_sorted = df_ref.sort_values("oscillator_strengths", ascending=False).head(top_k).reset_index(drop=True)
    
    # Ensure we only process valid peaks (non-NaN energies)
    df_ref_sorted = df_ref_sorted.dropna(subset=['exc_energies_ev'])
    
    n_targets = len(df_ref_sorted)
    
    if n_targets == 0:
        return [], []

    results_list = []

    for idx, df_pred in enumerate(dfs_pred):
        # Sort predicted peaks by energy for easier searching (optional but good practice)
        df_pred_sorted = df_pred.sort_values("exc_energies_ev").reset_index(drop=True)
        
        # Track which predicted indices have been "used" to prevent double counting
        used_pred_indices = set()
        
        matches_details = {
            "method_idx": idx,
            "matches": [],
            "energy_errors": [],
            "osc_errors": [],
            "unmatched_count": 0
        }
        
        total_error_energy = 0.0
        total_error_osc = 0.0
        
        # Iterate through the target (Reference) peaks
        for i, ref_row in df_ref_sorted.iterrows():
            ref_e = ref_row['exc_energies_ev']
            ref_f = ref_row['oscillator_strengths']
            
            best_dist = np.inf
            best_j = -1
            
            # Find the closest UNUSED predicted peak within tolerance
            for j, pred_row in df_pred_sorted.iterrows():
                if j in used_pred_indices:
                    continue
                
                pred_e = pred_row['exc_energies_ev']
                dist = abs(pred_e - ref_e)
                
                if dist < best_dist:
                    best_dist = dist
                    best_j = j
            
            # Check if the best match is within tolerance
            if best_dist <= tol_eV and best_j != -1:
                # Found a match!
                used_pred_indices.add(best_j)
                pred_row = df_pred_sorted.iloc[best_j]
                
                match_info = {
                    "ref_idx": i,
                    "ref_E": ref_e,
                    "ref_f": ref_f,
                    "pred_idx": best_j,
                    "pred_E": pred_row['exc_energies_ev'],
                    "pred_f": pred_row['oscillator_strengths'],
                    "delta_e": best_dist
                }
                matches_details["matches"].append(match_info)
                
                # Calculate Errors
                matches_details["energy_errors"].append(best_dist)
                matches_details["osc_errors"].append(abs(ref_f - pred_row['oscillator_strengths']))
                
                total_error_energy += best_dist
                total_error_osc += abs(ref_f - pred_row['oscillator_strengths'])
                
            else:
                # No match found within tolerance
                matches_details["unmatched_count"] += 1
                # Apply penalty for unmatched peaks
                matches_details["energy_errors"].append(unmatched_penalty)
                matches_details["osc_errors"].append(ref_f) # Or penalty? Usually just flag it.
                total_error_energy += unmatched_penalty
    
        # Calculate Metrics
        n_matched = len(matches_details["matches"])
        avg_err_E = np.mean(matches_details["energy_errors"]) if matches_details["energy_errors"] else float('nan')
        avg_err_f = np.mean(matches_details["osc_errors"]) if matches_details["osc_errors"] else float('nan')
        
        results = {
            "method_idx": idx,
            "n_targets": n_targets,
            "n_matched": n_matched,
            "match_ratio": n_matched / n_targets,
            "avg_energy_error_eV": avg_err_E,
            "avg_osc_error": avg_err_f,
            "unmatched_count": matches_details["unmatched_count"],
            "details": matches_details
        }
        
        results_list.append(results)
    
    return results_list, df_ref_sorted.to_dict('records')


def expand_gate(gate):
    """Expand one gate while preserving execution order."""
    if hasattr(gate, "gates"):
        return gate.gates
    return [gate]


def get_all_gates(circuit):
    expanded = []

    for op in circuit:
        if hasattr(op, "gates"):
            expanded.extend(op.gates)
        else:
            expanded.append(op)
    return expanded


def circuit_depth(circuit):
    # qubit -> last layer it was used in
    last_layer = defaultdict(int)

    depth = 0

    for op in circuit:
        ops = expand_gate(op)

        for g in ops:
            q = g.qubit_indices

            # compute when this gate can execute
            start = max((last_layer[i] for i in q), default=0)
            end = start + 1

            for i in q:
                last_layer[i] = end

            depth = max(depth, end)

    return depth


def get_fci_params(molecule_coords: list, basis: str, charge=0, multiplicity=1):
    """
    Calculate the number of electrons and spatial orbitals for a Full CI calculation
    given a molecule's coordinates and basis set.

    Parameters
    ----------
    molecule_coords : str or list
        Atom coordinates in PySCF format (e.g., "O 0 0 0; H 0 0 1.2" or list of tuples).
    basis : str
        Basis set name (e.g., 'sto-3g', '6-31g', 'cc-pvdz').
    charge : int, optional
        Total molecular charge (default 0).
    multiplicity : int, optional
        Spin multiplicity (2S+1) (default 1 for singlet).

    Returns
    -------
    dict
        A dictionary containing:
        - 'n_electrons': Total number of electrons.
        - 'n_orbitals': Total number of spatial orbitals.
        - 'n_spin_orbitals': Total number of spin orbitals (2 * n_orbitals).
        - 'charge': The input charge.
        - 'multiplicity': The input multiplicity.
        
    Raises
    ------
    RuntimeError
        If the SCF calculation fails to converge.
    """
    # 1. Define the Molecule
    mol = gto.M(
        atom=molecule_coords,
        basis=basis,
        charge=charge,
        spin=multiplicity - 1, # PySCF uses spin (2S) not multiplicity
        verbose=0 # Suppress output
    )

    # 2. Run Restricted Hartree-Fock (RHF) to get orbital count
    # Note: We run RHF just to get the basis set dimensions. 
    # If the system is open-shell, UHF might be needed, but orbital count is usually the same.
    mf = scf.RHF(mol).run()


    # 3. Extract parameters
    # mo_coeff shape is (n_basis, n_orbitals)
    n_spatial_orbitals = mf.mo_coeff.shape[1]
    n_electrons = sum(mol.nelec) # Returns tuple (alpha, beta), sum gives total
    n_spin_orbitals = n_spatial_orbitals * 2

    return {
        "n_electrons": n_electrons,
        "n_orbitals": n_spatial_orbitals,
        "n_spin_orbitals": n_spin_orbitals,
        "charge": charge,
        "multiplicity": multiplicity
    }


def count_singlet_states(n, m):
    """
    Calculate the number of singlet CSFs for a CAS(n, m) active space.
    
    Args:
        n (int): Number of active electrons (must be even).
        m (int): Number of active spatial orbitals.
        
    Returns:
        int: The number of singlet states.
        
    Raises:
        ValueError: If n is odd or n > 2*m.
    """
    if n < 0 or m < 0 or n > 2 * m or n % 2 != 0:
        raise ValueError("Invalid CAS(n, m) parameters for singlet calculation.")
    
    k = n // 2
    return (comb(m + 1, k) * comb(m + 1, k + 1)) // (m + 1)


def spectral_similarity(y1, y2, x):
    """
    Normalized spectral similarity (cosine analog):
    
        ∫ S1·S2 dE  /  √(∫ S1² dE · ∫ S2² dE)
    
    Returns value in [0, 1]:
        1.0 = identical spectra
        0.0 = completely disjoint
    """
    num = np.trapezoid(y1 * y2, x)
    den = np.sqrt(np.trapezoid(y1**2, x) * np.trapezoid(y2**2, x))
    if den == 0:
        return 0.0
    return num / den


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


def build_spectrum(energies, oscillator_strengths, x, sigma=0.2):
    """Build a Gaussian-broadened spectrum from discrete peaks."""
    y = np.zeros_like(x)
    for e, osc in zip(energies, oscillator_strengths):
        y += gaussian(x, e, osc, sigma)
    return y


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
    },
    "LiH": {
        "coords": [
            ("Li", 0.0, 0.0, 0.0),
            ("H", 0.0, 0.0, 1.595),
        ],
    },
    "HF": {
        "coords": [
            ("H", 0.0, 0.0, 0.0),
            ("F", 0.0, 0.0, 0.917),
        ],
    },
    "N2": {
        "coords": [
            ("N", 0.0, 0.0, -0.550),
            ("N", 0.0, 0.0, 0.550),
        ],
    },
    "CO": {
        "coords": [
            ("C", 0.0, 0.0, 0.0),
            ("O", 0.0, 0.0, 1.128),
        ],
    },
    "CH4": {
        "coords": [
            ("C", 0.0, 0.0, 0.0),
            ("H", 0.0, 0.0, 1.085),
            ("H", 0.0, 0.943, -0.362),
            ("H", 0.817, -0.471, -0.362),
            ("H", -0.817, -0.471, -0.362),
        ],
    },
    "NH3": {
        "coords": [
            ("N", 0.0, 0.0, 0.114),
            ("H", 0.0, 0.938, -0.342),
            ("H", 0.812, -0.469, -0.342),
            ("H", -0.812, -0.469, -0.342),
        ],
    },
    "BeH2": {
        "coords": [
            ("Be", 0.0, 0.0, 0.0),
            ("H", 0.0, 0.0, 1.330),
            ("H", 0.0, 0.0, -1.330),
        ],
    },
    "F2": {
        "coords": [
            ("F", 0.0, 0.0, -0.700),
            ("F", 0.0, 0.0, 0.700),
        ],
    },
    "R-methyloxirane": {
        "coords": [
            ("O",  0.8171, -0.7825, -0.2447),
            ("C", -1.5018,  0.1019, -0.1473),
            ("H", -1.3989,  0.3323, -1.2066),
            ("H", -2.0652, -0.8262, -0.0524),
            ("H", -2.0715,  0.8983,  0.3329),
            ("C", -0.1488, -0.0393,  0.4879),
            ("H", -0.1505, -0.2633,  1.5506),
            ("C",  1.0387,  0.6105, -0.0580),
            ("H",  0.9518,  1.2157, -0.9531),
            ("H",  1.8684,  0.8649,  0.5908),
        ],
    },
    "benzene": {
        "coords":[
            ("C", -1.210, 0.698, 0.004),
            ("C", -1.210, -0.698, 0.000),
            ("C", 0.000, -1.397, -0.003),
            ("C", 1.210, -0.698, -0.003),
            ("C", 1.210, 0.698, 0.000),
            ("C", 0.000, 1.397, -0.003),
            ("H", -2.164, 1.249, 0.007),
            ("H", -2.164, -1.249, 0.007),
            ("H", 0.000, -2.500, -0.005),
            ("H", 2.164, -1.249, -0.006),
            ("H", 2.164, 1.249, 0.000),
            ("H", 0.000, 2.500, 0.000),
        ]
    }
}


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
        scale_factor: float = None,
        # VQE Specific Inputs
        verbose: int = 0,
        casci_like: bool = False,
        calculator: Any = None,
        shots: int = None,
    ):
        """
        Initialize the workflow.
        
        Args:
            basis: Basis set string (e.g., 'aug-cc-pVDZ').
            molecule: molecule coordinates as list of tuples.
            num_active_orbitals: Number of active orbitals.
            num_active_electrons: Number of active electrons.
            num_states: Number of roots/states.
            verbose: Verbosity level.
            casci_like: Run DMDM like CASCI.
            calculator: qrunch calculator object for VQE
            shots: number of shots
        """
        self.casci_like = casci_like
        self.basis = basis
        self.molecule_list = molecule
        self.num_active_orbitals = num_active_orbitals
        self.num_active_electrons = num_active_electrons 
        self.num_states = num_states
        self.hartree_to_ev = 27.2114
        self.verbose = verbose
        self.calculator = calculator


        # Results Storage
        self._casci_results: Optional[Dict] = None
        self._vqe_results: Optional[Dict] = None
        self._dmdm_casci: Optional[Any] = None
        self._dmdm_vqe: Optional[Any] = None
                
        # Internal state flags
        self.scale_factor = scale_factor

        # Initial molecule setup
        if scale_factor is not None:
            self.molecule = scale_molecule(
                self.molecule_list,
                self.scale_factor,
                self.basis
            )
        else:
            self.molecule = gto.M(
                atom=molecule,
                basis=basis
            )

        self.shots = shots


    # ==========================
    # CLASSICAL (CASCI) PATH
    # ==========================

    def run_classical_casci_dmdm(self) -> Dict[str, Any]:
        """Run the classical CASCI workflow using PySCF."""

        start = time.time()
        process = psutil.Process(os.getpid())
        mem_before = process.memory_info().rss
        # 2. RHF
        mf = scf.RHF(self.molecule).run()

        # 3. CASCI
        mc = mcscf.CASCI(mf, ncas=self.num_active_orbitals, nelecas=self.num_active_electrons)    
        e_tot, e_cas, ci, mo, mo_energy = mc.kernel()
        self.casci_dmdm_ground_state = e_tot
        
        # 4. Integrals
        h_mo, _ = mc.get_h1eff()
        g_mo = mc.get_h2eff()
        g_mo = ao2mo.restore(1, g_mo, self.num_active_orbitals)

        rdm1, rdm2, rdm3, rdm4 = mc.fcisolver.make_rdm1234(ci, self.num_active_orbitals, self.num_active_electrons)

        # Dipole Integrals
        x_ao, y_ao, z_ao = self.molecule.intor('int1e_r', comp=3)
        cas_slice = slice(mc.ncore, mc.ncore + self.num_active_orbitals)
        C_cas = mo[:, cas_slice]
        
        x_cas = one_electron_integral_transform(C_cas, x_ao)
        y_cas = one_electron_integral_transform(C_cas, y_ao)
        z_cas = one_electron_integral_transform(C_cas, z_ao)
        MO_DM = [x_cas, y_cas, z_cas]

        mem_after = process.memory_info().rss
        self.mem_method = (mem_after - mem_before) / 1024 / 1024
        self.casci_dmdm_time_method = time.time() - start

        # Initialize DMDM
        dmdm = DMDM(
            h_mo,
            g_mo,
            0,
            self.num_active_orbitals,
            0,
            self.num_active_electrons,
            rdm1,
            rdm2=rdm2,
            rdm3=rdm3,
            rdm4=rdm4
        )

        # Get energies
        exc_energies = dmdm.get_excitation_energies() * self.hartree_to_ev

        # Get Oscillator Strengths
        osc_strengths = dmdm.get_oscillator_strength(MO_DM)

        mem_after2 = process.memory_info().rss
        self.mem_total = (mem_after2 - mem_before) / 1024 / 1024
        self.mem_dmdm = self.mem_method + (mem_after2 - mem_after) / 1024 / 1024
        self.casci_dmdm_time = time.time() - start
        self.casci_dmdm_dmdm = dmdm

        x, y, z = self.molecule.intor('int1e_cg_irxp', comp=3)

        MO_MG = [
            one_electron_integral_transform(C_cas, x),
            one_electron_integral_transform(C_cas, y),
            one_electron_integral_transform(C_cas, z)
        ]

        self.rotational_strengths = dmdm.get_rotational_strength(MO_DM, MO_MG)

        self._casci_dmdm_results = {
            'exc_energies_ev': exc_energies[:self.num_states],
            'oscillator_strengths': osc_strengths[:self.num_states],
            "rotational_strengths": self.rotational_strengths[:self.num_states]
        }
        if self.verbose > 0:
            print("Done with CASCI computations...\n")
        return self._casci_dmdm_results


    def run_classical_casci(self) -> Dict[str, Any]:
        """Run the classical CASCI workflow using PySCF."""

        start = time.time()
        process = psutil.Process(os.getpid())
        mem_before = process.memory_info().rss

        # 1. RHF
        mf = scf.RHF(self.molecule).run()

        # 2. CASCI Setup
        mc = mcscf.CASCI(mf, ncas=self.num_active_orbitals, nelecas=self.num_active_electrons)

        solver = fci.direct_spin1.FCI(self.molecule)
        solver = fci.addons.fix_spin_(solver, ss=0.0, shift=1.0)
        solver.spin = 0
        mc.fcisolver = solver
        mc.fcisolver.nroots = self.num_states

        # Run kernel
        e_tot, e_cas, ci, mo, mo_energy = mc.kernel()
        self.casci_ground_state = e_tot
        
        # 3. Integral Transformations
        cas_slice = slice(mc.ncore, mc.ncore + self.num_active_orbitals)
        C_cas = mo[:, cas_slice]
        
        # --- Electric Dipole Integrals (r) ---
        x_ao, y_ao, z_ao = self.molecule.intor('int1e_r', comp=3)
        x_cas = one_electron_integral_transform(C_cas, x_ao)
        y_cas = one_electron_integral_transform(C_cas, y_ao)
        z_cas = one_electron_integral_transform(C_cas, z_ao)
        MO_DM = [x_cas, y_cas, z_cas]
        
        # --- Magnetic Dipole Integrals (Angular Momentum r x p) ---
        mx_ao, my_ao, mz_ao = self.molecule.intor('int1e_cg_irxp', comp=3)
        mx_cas = one_electron_integral_transform(C_cas, mx_ao)
        my_cas = one_electron_integral_transform(C_cas, my_ao)
        mz_cas = one_electron_integral_transform(C_cas, mz_ao)
        MO_MG = [mx_cas, my_cas, mz_cas]

        # 4. Compute Excitation Energies, Oscillator Strengths, AND Rotational Strengths
        exc_energies_ev = []
        osc_strengths = []
        rot_strengths = []

        gs_ci = ci[0]
        gs_energy = e_cas[0] 

        for i in range(1, self.num_states):
            excited_ci = ci[i]
            excited_energy = e_cas[i]

            # Excitation Energy
            delta_e_hartree = excited_energy - gs_energy
            delta_e_ev = delta_e_hartree * self.hartree_to_ev
            exc_energies_ev.append(delta_e_ev)

            # Use trans_rdm1 ONLY (One-body operators only require 1-RDM)
            trans_rdm1 = mc.fcisolver.trans_rdm1(
                gs_ci, 
                excited_ci, 
                self.num_active_orbitals, 
                self.num_active_electrons
            )

            # --- Electric Transition Dipole Moment ---
            mu_x = np.einsum('pq,pq', trans_rdm1, x_cas)
            mu_y = np.einsum('pq,pq', trans_rdm1, y_cas)
            mu_z = np.einsum('pq,pq', trans_rdm1, z_cas)
            mu_sq = mu_x**2 + mu_y**2 + mu_z**2
            
            # Oscillator Strength
            f_val = (2.0 / 3.0) * delta_e_hartree * mu_sq
            osc_strengths.append(f_val)

            # --- Magnetic Transition Dipole Moment ---
            m_x = -np.einsum('pq,pq', trans_rdm1, mx_cas)
            m_y = -np.einsum('pq,pq', trans_rdm1, my_cas)
            m_z = -np.einsum('pq,pq', trans_rdm1, mz_cas)

            # Rotational Strength
            R_val = mu_x * m_x + mu_y * m_y + mu_z * m_z
            rot_strengths.append(R_val)
        
        mem_after = process.memory_info().rss
        self.mem_method = (mem_after - mem_before) / 1024 / 1024
        self.mem_total = (mem_after - mem_before) / 1024 / 1024
        self.casci_time = time.time() - start
        
        self._casci_results = {
            'exc_energies_ev': np.array(exc_energies_ev),
            'oscillator_strengths': np.array(osc_strengths),
            'rotational_strengths': np.array(rot_strengths),
        }
        
        if self.verbose > 0:
            print("Done with CASCI computations...\n")
            
        return self._casci_results

    # ==========================
    # QUANTUM (VQE) PATH
    # ==========================

    def run_quantum_vqe(self, use_noisy_backend: bool = False) -> Dict[str, Any]:
        """Run the VQE workflow using qrunch/qchem with optional noise simulation."""

        start = time.time()
        process = psutil.Process(os.getpid())
        mem_before = process.memory_info().rss
        # 1. Build Configuration (Unchanged)
        molecular_configuration = qc.build_molecular_configuration(molecule=self.molecule_list, basis_set=self.basis)
        
        num_inactive_orbs = molecular_configuration.number_of_alpha_electrons() - (self.num_active_electrons // 2)
        num_spatial_orbs = molecular_configuration.number_of_spatial_orbitals()
        num_virtual_orbs = num_spatial_orbs - num_inactive_orbs - self.num_active_orbitals

        print(f"  Inactive: {num_inactive_orbs}, Active: {self.num_active_orbitals}, Virtual: {num_virtual_orbs}")

        # Build Problem
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

        # 3. Run Calculation (Unchanged)
        result = self.calculator.calculate(ground_state_problem)
        self.vqe_ground_state = result.total_energy.value

        # 4. Extract Integrals (Unchanged)
        self.molecule = molecular_configuration.pyscf_molecule
        mo_coeffs = problem_builder._restricted_builder._calculate_molecular_orbitals(molecular_configuration).alpha.coefficients

        h_mo = one_electron_integral_transform(mo_coeffs, self.molecule.intor("int1e_kin") + self.molecule.intor("int1e_nuc"))
        g_mo = two_electron_integral_transform(mo_coeffs, self.molecule.intor("int2e"))
        effective_shots = self.shots if use_noisy_backend else (self.shots if self.shots else None)
        
        if use_noisy_backend and effective_shots is None:
            effective_shots = 4096 # Default to a reasonable shot count if user didn't specify
            print(f"  [WARN] Noisy backend requires shots. Defaulting to {effective_shots}.")

        
        if use_noisy_backend:
            print("  [INFO] Configuring Noisy Backend...")
            
            estimator_noisy = (
                qc.estimator_creator()
                .backend()
                .choose_backend()
                .local_qiskit_aer(method="density_matrix")
                .create()
            )

            estimator_exact = (
                qc.estimator_creator()
                .memory_restricted() # Default noiseless
                .with_precise_defaults()
                .create()
            )


            print(f"  [INFO] Noisy backend active. Shots: {self.shots}")
            rdm_calculator_noisy = ReducedDensityMatrixCalculator(estimator=estimator_noisy)
            rdm_calculator_exact = ReducedDensityMatrixCalculator(estimator=estimator_exact)

            rdm1 = rdm_calculator_noisy.calculate_1_rdm(circuit=result.final_circuit, shots=effective_shots)
            rdm2 = rdm_calculator_noisy.calculate_2_rdm(circuit=result.final_circuit, shots=effective_shots)
            # The noisy backend does not work with higher rdms
            rdm3 = rdm_calculator_exact.calculate_3_rdm(circuit=result.final_circuit, shots=None)
            rdm4 = rdm_calculator_exact.calculate_4_rdm(circuit=result.final_circuit, shots=None)
        else:
            print("  [INFO] Configuring Ideal (Noiseless) Backend...")
            estimator = (
                qc.estimator_creator()
                .memory_restricted()
                .with_precise_defaults()
                .create()
            )

            rdm_calculator = ReducedDensityMatrixCalculator(estimator=estimator)

            rdm1 = rdm_calculator.calculate_1_rdm(circuit=result.final_circuit, shots=None)
            rdm2 = rdm_calculator.calculate_2_rdm(circuit=result.final_circuit, shots=None)
            rdm3 = rdm_calculator.calculate_3_rdm(circuit=result.final_circuit, shots=None)
            rdm4 = rdm_calculator.calculate_4_rdm(circuit=result.final_circuit, shots=None)

        mem_after = process.memory_info().rss
        self.mem_method = (mem_after - mem_before) / 1024 / 1024
        self.vqe_time_method = time.time() - start
        mem_before2 = process.memory_info().rss
        # 6. DMDM Calculation (Unchanged logic, just using the new RDMs)
        start_dmdm = time.time()
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
            cas_slice = slice(num_inactive_orbs, num_inactive_orbs + self.num_active_orbitals)
            h_mo[cas_slice, cas_slice] = ground_state_problem.electronic_structure_integrals.one_body_core_hamiltonian.alpha_alpha
            g_mo[cas_slice, cas_slice, cas_slice, cas_slice] = ground_state_problem.electronic_structure_integrals.two_body_electron_repulsion_integrals.alpha_alpha

            # need to slice here since we get full space
            dmdm = DMDM(
                h_mo[cas_slice, cas_slice],
                g_mo[cas_slice, cas_slice, cas_slice, cas_slice],
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
        
        end_dmdm = time.time() - start_dmdm


        # 7. Dipole & Oscillator Strengths
        start_props = time.time()
        x, y, z = self.molecule.intor('int1e_r', comp=3)
        MO_DM = [
            one_electron_integral_transform(coefs, x),
            one_electron_integral_transform(coefs, y),
            one_electron_integral_transform(coefs, z)
        ]

        exc_energies_hartree = dmdm.get_excitation_energies()
        exc_energies_ev = exc_energies_hartree * self.hartree_to_ev
        osc_strengths = dmdm.get_oscillator_strength(MO_DM)

        mem_after2 = process.memory_info().rss
        self.mem_total = (mem_after2 - mem_before) / 1024 / 1024
        self.mem_dmdm = (mem_after2 - mem_before2) / 1024 / 1024
        self.vqe_time = time.time() - start
        self.vqe_dmdm = dmdm

        x, y, z = self.molecule.intor('int1e_cg_irxp', comp=3)

        MO_MG = [
            one_electron_integral_transform(coefs, x),
            one_electron_integral_transform(coefs, y),
            one_electron_integral_transform(coefs, z)
        ]

        rotational_strengths = dmdm.get_rotational_strength(MO_DM, MO_MG)
        self.rotational_strengths = rotational_strengths[exc_energies_ev > 1e-8]
        osc_strengths = osc_strengths[exc_energies_ev > 1e-8]
        exc_energies_ev = exc_energies_ev[exc_energies_ev > 1e-8]

        self._vqe_results = {
            'exc_energies_ev': np.append(exc_energies_ev[:self.num_states], [np.nan]*(self.num_states - exc_energies_ev[:self.num_states].shape[0])),
            'oscillator_strengths': np.append(osc_strengths[:self.num_states], [np.nan]*(self.num_states-osc_strengths[:self.num_states].shape[0])),
            "rotational_strengths": np.append(self.rotational_strengths[:self.num_states], [np.nan]*(self.num_states-self.rotational_strengths[:self.num_states].shape[0]))
        }
        
        if self.verbose > 0:
            print("Done with VQE computations...\n")
            
        return self._vqe_results


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