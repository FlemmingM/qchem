# Import PySCF, SlowQuant, DMDM
import pyscf
from pyscf import mcscf, mp, fci
import numpy as np

from dmdm.interface import DMDM

def one_electron_integral_transform(C: np.ndarray, int1e: np.ndarray) -> np.ndarray:
    return np.einsum("ai,bj,ab->ij", C, C, int1e, optimize=["einsum_path", (0, 2), (0, 1)])

def two_electron_integral_transform(C: np.ndarray, int2e: np.ndarray) -> np.ndarray:
    return np.einsum(
        "ai,bj,ck,dl,abcd->ijkl", C, C, C, C, int2e, optimize=["einsum_path", (0, 4), (0, 3), (0, 2), (0, 1)]
    )

AS_E, AS_O = 4, 4 # Active space: 4 electrons in 4 spatial orbitals
# Define molecule, basis set, and HF/MP2
basis = "sto-3g"
molecule_string = """
H 0.0 0.0 0.0
H 0.0 0.0 1.058354498
H 0.0 0.7937658735 0.0
H 0.0 0.7937658735 1.058354498
"""

molecule_string = """
O  0.000000   0.000000   0.000000
H  0.000000  -0.757160   0.586260
H  0.000000   0.757160   0.586260
"""


mol = pyscf.M(atom = molecule_string, basis=basis)
myhf = mol.RHF().run()
mymp = mp.MP2(myhf).run()
noons, natorbs = mcscf.addons.make_natural_orbitals(mymp)
cisolver = fci.FCI(mol, natorbs)
e, fcivec = cisolver.kernel()
print(f"E(FCI) = {e}")
rdm1, rdm2, rdm3, rdm4 = fci.rdm.make_dm1234('FCI4pdm_kern_sf', fcivec, fcivec, AS_O, AS_E)
rdm1, rdm2, rdm3, rdm4 = fci.rdm.reorder_dm1234(rdm1, rdm2, rdm3, rdm4)

h_mo = one_electron_integral_transform(natorbs, mol.intor("int1e_kin") + mol.intor("int1e_nuc"))
g_mo = two_electron_integral_transform(natorbs, mol.intor("int2e"))

dmdm = DMDM(
        h_mo,
        g_mo,
        0,
        AS_O,
        0,
        AS_E,
        rdm1,
        rdm2=rdm2,
        rdm3=rdm3,
        rdm4=rdm4,
        # print_level=3,
    )

x, y, z = mol.intor('int1e_r', comp=3)
MO_DM = [one_electron_integral_transform(natorbs, x), one_electron_integral_transform(natorbs, y), one_electron_integral_transform(natorbs, z)]

exc_energies = dmdm.get_excitation_energies()
osc_strengths = dmdm.get_oscillator_strength(MO_DM)

print("Excitation energies (a.u.): ", exc_energies)
print("Oscillator strengths: ", osc_strengths)
