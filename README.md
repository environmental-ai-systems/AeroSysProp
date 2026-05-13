# AeroSysProp

A Python library for predicting pure-component thermophysical properties of organic molecules relevant to atmospheric aerosol research. Given a SMILES string, AeroSysProp applies a range of published group-contribution methods to estimate boiling points, vapour pressures, critical properties, and liquid densities.

---

## Table of Contents

- [Background](#background)
- [Directory structure](#directory-structure)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Available methods](#available-methods)
- [Running the tests](#running-the-tests)
- [Docker](#docker)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Background

AeroSysProp builds on the earlier [UManSysProp](https://github.com/loftytopping/UManSysProp_pub) package. Its core purpose is to take a molecular structure (supplied as a SMILES string) and return property estimates using well-established group-contribution methods from the atmospheric chemistry literature. These estimates feed directly into aerosol partitioning models, droplet growth calculations, and related applications.

Molecular substructure matching is performed with [RDKit](https://www.rdkit.org/). Choice of predictive methods, and any modifications to them, is governed by the project team; the aim is always to provide a faithful implementation of the published technique, with any deviations clearly documented.

---

## Directory structure

```
AeroSysProp/
├── environment.yml                  # Conda environment definition
├── aerosysprop/
│   ├── __init__.py
│   ├── groups.py                    # SMARTS-based group counting and aggregation
│   ├── molar_volumes.py             # Molar volume helpers used by density methods
│   ├── data/                        # Group-contribution parameters and SMARTS patterns
│   │   ├── __init__.py              # Loads all .data, .smarts, and .matrix files at import
│   │   ├── *.smarts                 # SMARTS patterns for each method's functional groups
│   │   ├── *.data                   # Group contribution coefficients
│   │   └── *.matrix                 # Interaction matrices (e.g. Nannoolal)
│   ├── methods/
│   │   ├── boiling_points.py        # Joback & Reid, Stein & Brown, Nannoolal
│   │   ├── vapour_pressures.py      # Nannoolal, Myrdal & Yalkowsky, EVAPORATION, SIMPOL
│   │   ├── critical_properties.py   # Joback & Reid, Nannoolal
│   │   └── liquid_densities.py      # Girolami, Schroeder, Le Bas, Tyn & Calus
│   ├── examples/
│   │   ├── Pure_component_properties_example.py   # Worked example script
│   │   └── test.prop                              # Example SMILES input file
│   ├── tests/
│   │   ├── test_properties.py       # Regression test suite (run with pytest)
│   │   ├── fixtures.json            # Reference values for the regression tests
│   │   ├── generate_fixtures.py     # Script to regenerate fixtures.json
│   │   └── SMILES/
│   │       └── short_list.txt       # SMILES used for testing
│   └── Docker/
│       └── Dockerfile
```

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/loftytopping/AeroSysProp.git
cd AeroSysProp
```

### 2. Create the conda environment

The environment file installs Python, RDKit, NumPy, tqdm, and pytest:

```bash
conda env create -f environment.yml
conda activate AeroSysProp
```

### 3. Verify the installation

```bash
python -c "from aerosysprop.methods import boiling_points; print('OK')"
```

---

## Quick start

All methods accept an [RDKit](https://www.rdkit.org/) molecule object. The typical workflow is:

```python
from rdkit import Chem
from aerosysprop.methods import boiling_points, vapour_pressures, critical_properties, liquid_densities

# 1. Parse a SMILES string into an RDKit molecule
mol = Chem.MolFromSmiles('CC(O)OO')   # methyl hydroperoxide

# 2. Boiling point estimates (K)
bp_joback    = boiling_points.joback_and_reid(mol)
bp_stein     = boiling_points.stein_and_brown(mol)
bp_nannoolal = boiling_points.nannoolal(mol)

print(f'Joback & Reid:  {bp_joback:.2f} K')
print(f'Stein & Brown:  {bp_stein:.2f} K')
print(f'Nannoolal:      {bp_nannoolal:.2f} K')

# 3. Vapour pressure estimates (log10 atm) at 298.15 K
temperature = 298.15

vp_nannoolal = vapour_pressures.nannoolal(mol, temperature, bp_nannoolal)
vp_evap      = vapour_pressures.evaporation(mol, temperature)
vp_simpol    = vapour_pressures.simpol(mol, temperature)

print(f'Nannoolal VP:   {vp_nannoolal:.4f} log10(atm)')
print(f'EVAPORATION VP: {vp_evap:.4f} log10(atm)')
print(f'SIMPOL VP:      {vp_simpol:.4f} log10(atm)')

# 4. Critical properties (needed by some density methods)
cp = critical_properties.nannoolal(mol, bp_nannoolal)
print(f'Critical T:     {cp.temperature:.2f} K')
print(f'Critical P:     {cp.pressure:.4f} bar')

# 5. Liquid density (g/cc) at 298.15 K
density = liquid_densities.schroeder(mol, temperature, cp)
print(f'Schroeder rho:  {density.real:.4f} g/cc')
```

### Processing multiple compounds from a file

Input files are tab-separated with an index in column 1 and the SMILES string in column 2
(see [`aerosysprop/examples/test.prop`](aerosysprop/examples/test.prop)):

```
1    COO
2    CC(O)OO
3    CCCCC
```

```python
from rdkit import Chem
from aerosysprop.methods import boiling_points

with open('aerosysprop/examples/test.prop') as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 2:
            smiles = parts[1]
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                bp = boiling_points.nannoolal(mol)
                print(f'{smiles:40s}  Tb = {bp:.2f} K')
```

A fully worked example covering all property types is in
[`aerosysprop/examples/Pure_component_properties_example.py`](aerosysprop/examples/Pure_component_properties_example.py).

---

## Available methods

### Boiling points — `aerosysprop.methods.boiling_points`

| Function | Method | Returns |
|---|---|---|
| `joback_and_reid(mol)` | Joback & Reid (1987) | K |
| `stein_and_brown(mol)` | Stein & Brown (1994) | K |
| `nannoolal(mol)` | Nannoolal *et al.* (2004) | K |

### Vapour pressures — `aerosysprop.methods.vapour_pressures`

All functions return **log₁₀(P / atm)**. `T` is temperature in K.

| Function | Method | Requires boiling point? |
|---|---|---|
| `nannoolal(mol, T, bp)` | Nannoolal *et al.* (2008) | Yes |
| `myrdal_and_yalkowsky(mol, T, bp)` | Myrdal & Yalkowsky (1997) | Yes |
| `evaporation(mol, T)` | EVAPORATION (Compernolle *et al.*, 2011) | No |
| `evaporation2(mol, T)` | EVAPORATION with revised alcohol terms | No |
| `simpol(mol, T)` | SIMPOL.1 (Pankow & Asher, 2008) | No |

### Critical properties — `aerosysprop.methods.critical_properties`

Both functions return a named tuple `CriticalProperties(temperature, pressure, volume, compressibility)`.

| Function | Method |
|---|---|
| `joback_and_reid(mol, bp)` | Joback & Reid (1987) |
| `nannoolal(mol, bp)` | Nannoolal *et al.* (2004) |

### Liquid densities — `aerosysprop.methods.liquid_densities`

All functions return density in **g/cc**.

| Function | Method | Requires critical properties? |
|---|---|---|
| `girolami(mol)` | Girolami (1994) | No |
| `schroeder(mol, T, cp)` | Schroeder | Yes |
| `le_bas(mol, T, cp)` | Le Bas | Yes |
| `tyn_and_calus(mol, T, cp)` | Tyn & Calus | Yes |

---

## Running the tests

The regression test suite checks all methods against pre-computed reference values across 21 representative atmospheric compounds:

```bash
conda activate AeroSysProp
pytest aerosysprop/tests/test_properties.py -v
```

Expected output: **399 passed**.

### Regenerating reference values

If you intentionally change a method's implementation, regenerate `fixtures.json` and commit it alongside your code change:

```bash
python aerosysprop/tests/generate_fixtures.py
```

---

## Docker

A Dockerfile is provided for containerised use:

```bash
# Build (requires SSH access to the GitHub repo during the build)
docker build --ssh github=~/.ssh/id_rsa_github -t aerosysprop aerosysprop/Docker/

# Run interactively
docker run -it aerosysprop bash
```

Inside the container, activate the environment and use the library as normal:

```bash
conda activate AeroSysProp
python -c "from aerosysprop.methods import boiling_points; print('ready')"
```

---

## Contributing

Contributions are welcome. Please open an issue before submitting a pull request so that proposed method changes can be discussed with the project team. Any modification to an existing method should be accompanied by:

1. A clear description of what differs from the published technique.
2. Updated reference values (`fixtures.json`) generated by running `generate_fixtures.py`.
3. Passing tests (`pytest aerosysprop/tests/test_properties.py`).

---

## Citation

If you use AeroSysProp in published work, please cite the repository and the original papers for the methods you use (see the tables above).

---

## License

AeroSysProp is distributed under the GNU General Public License v3. See [LICENSE](LICENSE) for details.

---

## Contact

David Topping — david.topping@manchester.ac.uk  
University of Manchester
