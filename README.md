# Documentation

## Installation

### 1. Create a micromamba environment

```bash
micromamba create -n qchem python=3.11 -y
micromamba activate qchem
```

### 2. Install Python dependencies
``` 
pip install pyscf numpy pandas matplotlib psutil
```

### 3. Install the DMDM package from GitHub
```
git clone https://github.com/HQC2/DensityMatrixDrivenModule.git
cd DensityMatrixDrivenModule
pip install .
```

### 4. Install qrunch (for VQE support)
You will need a licence from kvantify.com.

```
pip install qrunch --extra-index-url $KVANTIFY_INDEX_URL
```


## Quick Start
```
from dmdm_workflow import DMDMWorkflow

# Example: Water molecule
molecule = [
    ("O", 0.0, 0.0, 0.1035),
    ("H", 0.0, 0.7956, -0.4640),
    ("H", 0.0, -0.7956, -0.4640),
]
```

## Run Classical CASCI Workflow
```
workflow = DMDMWorkflow(
    basis='sto-3g',
    molecule=molecule,
    num_active_orbitals=4,      # Number of active orbitals
    num_active_electrons=4,     # Number of active electrons
    num_states=20,              # How many excited states to compute
    verbose=1                   # Set > 0 for progress output
)

results = workflow.run_classical_casci()
print(results['exc_energies_ev'])       # Excitation energies in eV
print(results['oscillator_strengths'])  # Oscillator strengths
print(results['rotational_strengths'])  # Rotational strengths
```

## Run Quantum VQE Workflow
```
calculator = (
    qc.calculator_creator()
        .vqe()
        .iterative()
        .standard()
        .with_options(
                options=qc.options.IterativeVqeOptions(max_iterations=500)
        )
        .create()
)

workflow = DMDMWorkflow(
    basis='sto-3g',
    molecule=molecule,
    num_active_orbitals=6,
    num_active_electrons=8,
    num_states=5,
    calculator=calculator,      # Required for VQE path
    shots=4096,                 # Shot count for noisy simulations
    verbose=1
)
```

## Noiseless simulation
```
results = workflow.run_quantum_vqe(use_noisy_backend=False)
```

## Noisy simulation (requires shots parameter)
```
results = workflow.run_quantum_vqe(use_noisy_backend=True)
```