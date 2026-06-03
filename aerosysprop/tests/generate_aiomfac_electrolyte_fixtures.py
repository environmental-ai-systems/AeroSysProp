"""Generate Fortran AIOMFAC reference outputs for Phase 2 (electrolyte)
regression tests.

For each curated mixture (neutrals + salts), this script:
  1. Decomposes each neutral SMILES via groups.aiomfac.
  2. Looks up each salt in the AIOMFAC catalog.
  3. Builds an AIOMFAC Fortran input with mass fractions.
  4. Runs the Fortran executable.
  5. Parses water/organic gamma_x AND ion gamma_m from the output.
  6. Aggregates into aiomfac_electrolyte_fixtures.json.

Mole-fraction-basis gamma reported by Fortran for neutrals (a_coeff_x);
molality-basis gamma for ions (a_coeff_m).
"""

import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from aerosysprop import data, groups
from rdkit import Chem


HERE = Path(__file__).parent
FIXTURE_FILE = HERE / 'aiomfac_electrolyte_fixtures.json'

FORTRAN_DIR = Path('/Users/davidtopping/Documents/GitHub/AIOMFAC/FortranCode')
FORTRAN_EXE = FORTRAN_DIR / 'AIOMFAC-web.out'

TEMPERATURE_K = 298.15
TOPSUBNO = 265
EXTENSION_IDS = {281, 282, 283, 284, 285}


# Each case: (case_id, [(neutral_name, smiles), ...], [salt_name, ...],
#             [{component_name: mass_fraction}, ...])
# Water is mandatory as the first neutral (Fortran convention).
CASES = [
    (
        'water_NaCl', [('water', 'O')], ['NaCl'],
        [{'NaCl': w} for w in (0.01, 0.05, 0.10, 0.15, 0.20, 0.25)],
    ),
    (
        'water_KCl', [('water', 'O')], ['KCl'],
        [{'KCl': w} for w in (0.01, 0.05, 0.10, 0.15, 0.20, 0.25)],
    ),
    (
        'water_NH4Cl', [('water', 'O')], ['NH4Cl'],
        [{'NH4Cl': w} for w in (0.01, 0.05, 0.10, 0.15, 0.20, 0.25)],
    ),
    (
        'water_NaNO3', [('water', 'O')], ['NaNO3'],
        [{'NaNO3': w} for w in (0.01, 0.05, 0.10, 0.20, 0.30, 0.40)],
    ),
    (
        'water_NH4NO3', [('water', 'O')], ['NH4NO3'],
        [{'NH4NO3': w} for w in (0.01, 0.05, 0.10, 0.20, 0.30, 0.40)],
    ),
    (
        'water_Na2SO4', [('water', 'O')], ['Na2SO4'],
        [{'Na2SO4': w} for w in (0.01, 0.05, 0.10, 0.15, 0.20)],
    ),
    (
        'water_NH42SO4', [('water', 'O')], ['NH42SO4'],  # (NH4)2SO4
        [{'NH42SO4': w} for w in (0.01, 0.05, 0.10, 0.20, 0.30)],
    ),
    (
        'water_HCl', [('water', 'O')], ['HCl'],
        [{'HCl': w} for w in (0.01, 0.03, 0.05, 0.10, 0.15)],
    ),
    (
        'water_HNO3', [('water', 'O')], ['HNO3'],
        [{'HNO3': w} for w in (0.01, 0.03, 0.05, 0.10, 0.20)],
    ),
    # H2SO4 requires the bisulfate dissociation equilibrium (HSO4- <-> H+ + SO4--)
    # which Fortran solves via SubModDissociationEquil.f90. The current Python
    # port treats H2SO4 as fully dissociated to 2 H+ + SO4-- so it would not
    # match Fortran's split between HSO4- and SO4--. Mark as Phase 3.
    # (excluded: 'water_H2SO4')
    # ---- water + organic + salt ----
    (
        'water_ethanol_NaCl',
        [('water', 'O'), ('ethanol', 'CCO')], ['NaCl'],
        [
            {'ethanol': 0.10, 'NaCl': 0.02},
            {'ethanol': 0.20, 'NaCl': 0.05},
            {'ethanol': 0.30, 'NaCl': 0.05},
            {'ethanol': 0.10, 'NaCl': 0.10},
            {'ethanol': 0.50, 'NaCl': 0.02},
        ],
    ),
    (
        'water_glucose_NaCl',
        [('water', 'O'), ('glucose', 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O')], ['NaCl'],
        [
            {'glucose': 0.05, 'NaCl': 0.02},
            {'glucose': 0.10, 'NaCl': 0.05},
            {'glucose': 0.20, 'NaCl': 0.05},
        ],
    ),
    (
        'water_ethanol_NH42SO4',
        [('water', 'O'), ('ethanol', 'CCO')], ['NH42SO4'],
        [
            {'ethanol': 0.10, 'NH42SO4': 0.05},
            {'ethanol': 0.20, 'NH42SO4': 0.10},
            {'ethanol': 0.30, 'NH42SO4': 0.05},
        ],
    ),
    (
        'water_aceticacid_NaCl',
        [('water', 'O'), ('acetic_acid', 'CC(=O)O')], ['NaCl'],
        [
            {'acetic_acid': 0.10, 'NaCl': 0.05},
            {'acetic_acid': 0.30, 'NaCl': 0.05},
            {'acetic_acid': 0.50, 'NaCl': 0.05},
        ],
    ),
    (
        'water_succinicacid_NH42SO4',
        [('water', 'O'), ('succinic_acid', 'OC(=O)CCC(=O)O')], ['NH42SO4'],
        [
            {'succinic_acid': 0.05, 'NH42SO4': 0.05},
            {'succinic_acid': 0.10, 'NH42SO4': 0.10},
        ],
    ),
    # ---- atmospheric SOA proxy: organic + sulfate aerosol ----
    (
        'water_glycerol_NaCl',   # Triol + chloride: glycerol now correctly decomposes
        [('water', 'O'), ('glycerol', 'OCC(O)CO')], ['NaCl'],
        [
            {'glycerol': 0.10, 'NaCl': 0.03},
            {'glycerol': 0.20, 'NaCl': 0.05},
            {'glycerol': 0.30, 'NaCl': 0.05},
            {'glycerol': 0.10, 'NaCl': 0.10},
        ],
    ),
    (
        'water_pinonic_NH42SO4',
        [('water', 'O'), ('pinonic_acid', 'CC(=O)C1CCC(C1(C)C)CC(=O)O')], ['NH42SO4'],
        [
            {'pinonic_acid': 0.05, 'NH42SO4': 0.10},
            {'pinonic_acid': 0.10, 'NH42SO4': 0.15},
            {'pinonic_acid': 0.15, 'NH42SO4': 0.05},
        ],
    ),
    # ---- mixed-salt cases ----
    (
        'water_NaCl_KCl',
        [('water', 'O')], ['NaCl', 'KCl'],
        [
            {'NaCl': 0.03, 'KCl': 0.03},
            {'NaCl': 0.05, 'KCl': 0.10},
            {'NaCl': 0.10, 'KCl': 0.05},
        ],
    ),
    # ---- multi-salt: cross-coupled cation/anion pairs -----------------------
    (
        'water_NaCl_Na2SO4',  # Na+ with both Cl- and SO4--
        [('water', 'O')], ['NaCl', 'Na2SO4'],
        [
            {'NaCl': 0.03, 'Na2SO4': 0.02},
            {'NaCl': 0.05, 'Na2SO4': 0.05},
            {'NaCl': 0.10, 'Na2SO4': 0.05},
            {'NaCl': 0.02, 'Na2SO4': 0.10},
        ],
    ),
    (
        'water_NH42SO4_NH4NO3',   # canonical atmospheric ammonium sulfate-nitrate
        [('water', 'O')], ['NH42SO4', 'NH4NO3'],
        [
            {'NH42SO4': 0.05, 'NH4NO3': 0.05},
            {'NH42SO4': 0.10, 'NH4NO3': 0.10},
            {'NH42SO4': 0.20, 'NH4NO3': 0.05},
            {'NH42SO4': 0.05, 'NH4NO3': 0.20},
        ],
    ),
    (
        'water_NaCl_NH42SO4',   # 2 cations x 2 anions, all four pairs activated
        [('water', 'O')], ['NaCl', 'NH42SO4'],
        [
            {'NaCl': 0.03, 'NH42SO4': 0.03},
            {'NaCl': 0.05, 'NH42SO4': 0.10},
            {'NaCl': 0.10, 'NH42SO4': 0.05},
        ],
    ),
    (
        'water_KNO3_NaCl',
        [('water', 'O')], ['KNO3', 'NaCl'],
        [
            {'KNO3': 0.05, 'NaCl': 0.03},
            {'KNO3': 0.10, 'NaCl': 0.05},
            {'KNO3': 0.05, 'NaCl': 0.10},
        ],
    ),
    (
        'water_KCl_Na2SO4',
        [('water', 'O')], ['KCl', 'Na2SO4'],
        [
            {'KCl': 0.05, 'Na2SO4': 0.05},
            {'KCl': 0.10, 'Na2SO4': 0.10},
        ],
    ),
    # ---- multi-salt: 3+ salts -----------------------------------------------
    (
        'water_NaCl_KCl_NaNO3',
        [('water', 'O')], ['NaCl', 'KCl', 'NaNO3'],
        [
            {'NaCl': 0.03, 'KCl': 0.03, 'NaNO3': 0.03},
            {'NaCl': 0.05, 'KCl': 0.05, 'NaNO3': 0.05},
            {'NaCl': 0.02, 'KCl': 0.08, 'NaNO3': 0.05},
        ],
    ),
    (
        'water_NH42SO4_NH4Cl_NH4NO3',  # all-ammonium with three anions
        [('water', 'O')], ['NH42SO4', 'NH4Cl', 'NH4NO3'],
        [
            {'NH42SO4': 0.05, 'NH4Cl': 0.03, 'NH4NO3': 0.05},
            {'NH42SO4': 0.10, 'NH4Cl': 0.05, 'NH4NO3': 0.10},
            {'NH42SO4': 0.15, 'NH4Cl': 0.02, 'NH4NO3': 0.03},
        ],
    ),
    (
        'water_NaCl_NH42SO4_KCl_NaNO3',   # 4 salts, 3 cations, 3 anions
        [('water', 'O')], ['NaCl', 'NH42SO4', 'KCl', 'NaNO3'],
        [
            {'NaCl': 0.03, 'NH42SO4': 0.05, 'KCl': 0.02, 'NaNO3': 0.03},
            {'NaCl': 0.05, 'NH42SO4': 0.10, 'KCl': 0.03, 'NaNO3': 0.02},
        ],
    ),
    # ---- divalent cations (Mg2+, Ca2+) --------------------------------------
    (
        'water_MgCl2', [('water', 'O')], ['MgCl2'],
        [{'MgCl2': w} for w in (0.01, 0.05, 0.10, 0.20)],
    ),
    (
        'water_Ca(Cl)2', [('water', 'O')], ['Ca(Cl)2'],
        [{'Ca(Cl)2': w} for w in (0.01, 0.05, 0.10, 0.20)],
    ),
    (
        'water_MgSO4', [('water', 'O')], ['MgSO4'],
        [{'MgSO4': w} for w in (0.01, 0.05, 0.10, 0.15)],
    ),
    (
        'water_MgCl2_Ca(Cl)2', [('water', 'O')], ['MgCl2', 'Ca(Cl)2'],
        [
            {'MgCl2': 0.03, 'Ca(Cl)2': 0.03},
            {'MgCl2': 0.05, 'Ca(Cl)2': 0.10},
        ],
    ),
    (
        'water_NH42SO4_MgSO4',
        [('water', 'O')], ['NH42SO4', 'MgSO4'],
        [
            {'NH42SO4': 0.05, 'MgSO4': 0.03},
            {'NH42SO4': 0.10, 'MgSO4': 0.05},
        ],
    ),
    # ---- multi-organic + single salt ----------------------------------------
    (
        'water_ethanol_butanediol_NaCl',
        [('water', 'O'), ('ethanol', 'CCO'), ('1,4-butanediol', 'OCCCCO')],
        ['NaCl'],
        [
            {'ethanol': 0.10, '1,4-butanediol':0.10, 'NaCl': 0.05},
            {'ethanol': 0.20, '1,4-butanediol':0.10, 'NaCl': 0.05},
            {'ethanol': 0.10, '1,4-butanediol':0.20, 'NaCl': 0.03},
        ],
    ),
    (
        'water_ethanol_aceticacid_NH42SO4',
        [('water', 'O'), ('ethanol', 'CCO'), ('acetic_acid', 'CC(=O)O')],
        ['NH42SO4'],
        [
            {'ethanol': 0.10, 'acetic_acid': 0.10, 'NH42SO4': 0.05},
            {'ethanol': 0.20, 'acetic_acid': 0.10, 'NH42SO4': 0.10},
            {'ethanol': 0.30, 'acetic_acid': 0.20, 'NH42SO4': 0.05},
        ],
    ),
    (
        'water_glucose_aceticacid_NH42SO4',
        [('water', 'O'),
         ('glucose', 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'),
         ('acetic_acid', 'CC(=O)O')],
        ['NH42SO4'],
        [
            {'glucose': 0.05, 'acetic_acid': 0.05, 'NH42SO4': 0.05},
            {'glucose': 0.10, 'acetic_acid': 0.10, 'NH42SO4': 0.10},
            {'glucose': 0.15, 'acetic_acid': 0.05, 'NH42SO4': 0.05},
        ],
    ),
    (
        'water_4organic_NaCl',  # 4 organics + salt
        [('water', 'O'), ('ethanol', 'CCO'), ('1-propanol', 'CCCO'),
         ('ethylene_glycol', 'OCCO')],
        ['NaCl'],
        [
            {'ethanol': 0.10, '1-propanol': 0.05, 'ethylene_glycol': 0.10, 'NaCl': 0.03},
            {'ethanol': 0.15, '1-propanol': 0.10, 'ethylene_glycol': 0.15, 'NaCl': 0.05},
        ],
    ),
    # ---- realistic atmospheric SOA + multi-salt -----------------------------
    (
        'urban_aerosol_proxy',   # ethanol + acetic + ammonium sulfate-nitrate
        [('water', 'O'), ('ethanol', 'CCO'), ('acetic_acid', 'CC(=O)O')],
        ['NH42SO4', 'NH4NO3'],
        [
            {'ethanol': 0.10, 'acetic_acid': 0.10, 'NH42SO4': 0.05, 'NH4NO3': 0.05},
            {'ethanol': 0.15, 'acetic_acid': 0.15, 'NH42SO4': 0.10, 'NH4NO3': 0.05},
            {'ethanol': 0.05, 'acetic_acid': 0.10, 'NH42SO4': 0.10, 'NH4NO3': 0.10},
        ],
    ),
    (
        'marine_SOA_seasalt_sulfate',   # alpha-pinene SOA + sugar + sea-salt + sulfate
        [('water', 'O'),
         ('pinonic_acid', 'CC(=O)C1CCC(C1(C)C)CC(=O)O'),
         ('glucose', 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O')],
        ['NH42SO4', 'NaCl'],
        [
            {'pinonic_acid': 0.05, 'glucose': 0.05, 'NH42SO4': 0.05, 'NaCl': 0.03},
            {'pinonic_acid': 0.10, 'glucose': 0.05, 'NH42SO4': 0.10, 'NaCl': 0.02},
            {'pinonic_acid': 0.05, 'glucose': 0.10, 'NH42SO4': 0.05, 'NaCl': 0.05},
        ],
    ),
    (
        'water_succinic_butanediol_Na2SO4_KCl',
        [('water', 'O'),
         ('succinic_acid', 'OC(=O)CCC(=O)O'),
         ('1,4-butanediol', 'OCCCCO')],
        ['Na2SO4', 'KCl'],
        [
            {'succinic_acid': 0.05, '1,4-butanediol':0.05, 'Na2SO4': 0.05, 'KCl': 0.03},
            {'succinic_acid': 0.10, '1,4-butanediol':0.10, 'Na2SO4': 0.05, 'KCl': 0.05},
        ],
    ),
    (
        'biomass_burning_proxy',   # levoglucosan + acetic acid + sea-salt + sulfate
        [('water', 'O'),
         ('levoglucosan', 'OC[C@H]1OC2OC1C(O)C2O'),
         ('acetic_acid', 'CC(=O)O')],
        ['NaCl', 'NH42SO4'],
        [
            {'levoglucosan': 0.05, 'acetic_acid': 0.05, 'NaCl': 0.03, 'NH42SO4': 0.05},
            {'levoglucosan': 0.10, 'acetic_acid': 0.10, 'NaCl': 0.02, 'NH42SO4': 0.05},
            {'levoglucosan': 0.15, 'acetic_acid': 0.05, 'NaCl': 0.05, 'NH42SO4': 0.10},
        ],
    ),
    (
        'alpha_pinene_SOA_mixed_inorg',   # pinic + pinonic + sulfate + nitrate
        [('water', 'O'),
         ('pinic_acid', 'CC1(C)C(C(=O)O)CCC1CC(=O)O'),
         ('pinonic_acid', 'CC(=O)C1CCC(C1(C)C)CC(=O)O')],
        ['NH42SO4', 'NH4NO3'],
        [
            {'pinic_acid': 0.05, 'pinonic_acid': 0.05, 'NH42SO4': 0.05, 'NH4NO3': 0.05},
            {'pinic_acid': 0.10, 'pinonic_acid': 0.05, 'NH42SO4': 0.10, 'NH4NO3': 0.05},
            {'pinic_acid': 0.05, 'pinonic_acid': 0.10, 'NH42SO4': 0.05, 'NH4NO3': 0.10},
        ],
    ),
    # ---- edge concentrations ------------------------------------------------
    (
        'water_NaCl_very_dilute',  # 0.001 to 0.01 mass-fraction (near infinite dilution)
        [('water', 'O')], ['NaCl'],
        [{'NaCl': w} for w in (0.0001, 0.0010, 0.0050)],
    ),
    (
        'water_NH42SO4_high_conc',  # high concentration up to ~40%
        [('water', 'O')], ['NH42SO4'],
        [{'NH42SO4': w} for w in (0.25, 0.30, 0.35, 0.40)],
    ),
    (
        'water_NaCl_at_saturation',
        [('water', 'O')], ['NaCl'],
        [{'NaCl': w} for w in (0.20, 0.23, 0.25)],
    ),
    (
        'organic_rich_low_water',   # low water, organic-rich + small salt
        [('water', 'O'), ('ethanol', 'CCO')], ['NaCl'],
        [
            {'ethanol': 0.50, 'NaCl': 0.02},
            {'ethanol': 0.70, 'NaCl': 0.02},
            {'ethanol': 0.80, 'NaCl': 0.01},
        ],
    ),
    (
        'water_dilute_multi_salt',   # dilute mix of multiple salts
        [('water', 'O')], ['NaCl', 'KCl', 'NaNO3'],
        [
            {'NaCl': 0.001, 'KCl': 0.001, 'NaNO3': 0.001},
            {'NaCl': 0.005, 'KCl': 0.005, 'NaNO3': 0.005},
        ],
    ),
    # ====================================================================
    # High-organic-count stress tests (5, 10, 20 organics + 1-2 salts).
    # Each organic mass fraction is intentionally modest so we stay in a
    # realistic atmospheric-aerosol concentration regime.
    # ====================================================================
    (
        '5_organic_NaCl',
        [('water', 'O'),
         ('ethanol',         'CCO'),
         ('glycerol',        'OCC(O)CO'),
         ('acetic_acid',     'CC(=O)O'),
         ('glucose',         'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'),
         ('levoglucosan',    'OC[C@H]1OC2OC1C(O)C2O')],
        ['NaCl'],
        [
            {'ethanol': 0.05, 'glycerol': 0.05, 'acetic_acid': 0.05,
             'glucose': 0.05, 'levoglucosan': 0.05, 'NaCl': 0.03},
            {'ethanol': 0.10, 'glycerol': 0.05, 'acetic_acid': 0.05,
             'glucose': 0.03, 'levoglucosan': 0.02, 'NaCl': 0.05},
        ],
    ),
    (
        '5_organic_NH42SO4_NH4NO3',
        [('water', 'O'),
         ('methanol',        'CO'),
         ('ethanol',         'CCO'),
         ('1-propanol',      'CCCO'),
         ('ethylene_glycol', 'OCCO'),
         ('acetic_acid',     'CC(=O)O')],
        ['NH42SO4', 'NH4NO3'],
        [
            {'methanol': 0.03, 'ethanol': 0.03, '1-propanol': 0.03,
             'ethylene_glycol': 0.05, 'acetic_acid': 0.05,
             'NH42SO4': 0.05, 'NH4NO3': 0.05},
            {'methanol': 0.05, 'ethanol': 0.05, '1-propanol': 0.05,
             'ethylene_glycol': 0.10, 'acetic_acid': 0.05,
             'NH42SO4': 0.10, 'NH4NO3': 0.05},
        ],
    ),
    (
        '10_organic_NaCl',  # alcohols + polyols + sugars + salt
        [('water', 'O'),
         ('methanol',        'CO'),
         ('ethanol',         'CCO'),
         ('1-propanol',      'CCCO'),
         ('1-butanol',       'CCCCO'),
         ('ethylene_glycol', 'OCCO'),
         ('1,4-butanediol',  'OCCCCO'),
         ('glycerol',        'OCC(O)CO'),
         ('erythritol',      'OCC(O)C(O)CO'),
         ('xylitol',         'OCC(O)C(O)C(O)CO'),
         ('glucose',         'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O')],
        ['NaCl'],
        [
            {'methanol': 0.02, 'ethanol': 0.02, '1-propanol': 0.02,
             '1-butanol': 0.02, 'ethylene_glycol': 0.03,
             '1,4-butanediol': 0.03, 'glycerol': 0.03,
             'erythritol': 0.03, 'xylitol': 0.03, 'glucose': 0.03,
             'NaCl': 0.03},
            {'methanol': 0.03, 'ethanol': 0.03, '1-propanol': 0.03,
             '1-butanol': 0.02, 'ethylene_glycol': 0.05,
             '1,4-butanediol': 0.03, 'glycerol': 0.05,
             'erythritol': 0.02, 'xylitol': 0.02, 'glucose': 0.05,
             'NaCl': 0.05},
        ],
    ),
    (
        '10_organic_atm_realistic',   # SOA + biomass + acids + ammonium sulfate
        [('water', 'O'),
         ('pinic_acid',      'CC1(C)C(C(=O)O)CCC1CC(=O)O'),
         ('pinonic_acid',    'CC(=O)C1CCC(C1(C)C)CC(=O)O'),
         ('glucose',         'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'),
         ('levoglucosan',    'OC[C@H]1OC2OC1C(O)C2O'),
         ('acetic_acid',     'CC(=O)O'),
         ('ethanol',         'CCO'),
         ('glycerol',        'OCC(O)CO'),
         ('succinic_acid',   'OC(=O)CCC(=O)O'),
         ('glutaric_acid',   'OC(=O)CCCC(=O)O'),
         ('malic_acid',      'OC(=O)CC(O)C(=O)O')],
        ['NH42SO4'],
        [
            {'pinic_acid': 0.02, 'pinonic_acid': 0.02, 'glucose': 0.03,
             'levoglucosan': 0.02, 'acetic_acid': 0.03, 'ethanol': 0.03,
             'glycerol': 0.03, 'succinic_acid': 0.02, 'glutaric_acid': 0.02,
             'malic_acid': 0.02, 'NH42SO4': 0.05},
            {'pinic_acid': 0.03, 'pinonic_acid': 0.03, 'glucose': 0.05,
             'levoglucosan': 0.03, 'acetic_acid': 0.05, 'ethanol': 0.05,
             'glycerol': 0.05, 'succinic_acid': 0.03, 'glutaric_acid': 0.03,
             'malic_acid': 0.03, 'NH42SO4': 0.10},
        ],
    ),
    (
        '20_organic_NaCl',  # broad organic survey + chloride salt
        [('water', 'O'),
         ('methanol',         'CO'),
         ('ethanol',          'CCO'),
         ('1-propanol',       'CCCO'),
         ('1-butanol',        'CCCCO'),
         ('glycerol',         'OCC(O)CO'),
         ('ethylene_glycol',  'OCCO'),
         ('glucose',          'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'),
         ('levoglucosan',     'OC[C@H]1OC2OC1C(O)C2O'),
         ('sucrose',          'OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O'),
         ('formic_acid',      'OC=O'),
         ('acetic_acid',      'CC(=O)O'),
         ('propionic_acid',   'CCC(=O)O'),
         ('oxalic_acid',      'OC(=O)C(=O)O'),
         ('malonic_acid',     'OC(=O)CC(=O)O'),
         ('succinic_acid',    'OC(=O)CCC(=O)O'),
         ('glutaric_acid',    'OC(=O)CCCC(=O)O'),
         ('adipic_acid',      'OC(=O)CCCCC(=O)O'),
         ('malic_acid',       'OC(=O)CC(O)C(=O)O'),
         ('citric_acid',      'OC(=O)CC(O)(C(=O)O)CC(=O)O'),
         ('lactic_acid',      'CC(O)C(=O)O')],
        ['NaCl'],
        [
            # 20 organics × 0.01-0.02 each plus salt; water ~ 0.6
            dict(zip(
                ['methanol','ethanol','1-propanol','1-butanol','glycerol',
                 'ethylene_glycol','glucose','levoglucosan','sucrose',
                 'formic_acid','acetic_acid','propionic_acid',
                 'oxalic_acid','malonic_acid','succinic_acid',
                 'glutaric_acid','adipic_acid','malic_acid',
                 'citric_acid','lactic_acid','NaCl'],
                [0.01]*20 + [0.05],
            )),
            dict(zip(
                ['methanol','ethanol','1-propanol','1-butanol','glycerol',
                 'ethylene_glycol','glucose','levoglucosan','sucrose',
                 'formic_acid','acetic_acid','propionic_acid',
                 'oxalic_acid','malonic_acid','succinic_acid',
                 'glutaric_acid','adipic_acid','malic_acid',
                 'citric_acid','lactic_acid','NaCl'],
                [0.02]*20 + [0.03],
            )),
        ],
    ),
    (
        '20_organic_atm_SOA_sulfate_nitrate',
        # Realistic atmospheric SOA + secondary inorganic mix
        [('water', 'O'),
         ('pinic_acid',       'CC1(C)C(C(=O)O)CCC1CC(=O)O'),
         ('pinonic_acid',     'CC(=O)C1CCC(C1(C)C)CC(=O)O'),
         ('succinic_acid',    'OC(=O)CCC(=O)O'),
         ('glutaric_acid',    'OC(=O)CCCC(=O)O'),
         ('adipic_acid',      'OC(=O)CCCCC(=O)O'),
         ('glycerol',         'OCC(O)CO'),
         ('glucose',          'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'),
         ('levoglucosan',     'OC[C@H]1OC2OC1C(O)C2O'),
         ('acetic_acid',      'CC(=O)O'),
         ('lactic_acid',      'CC(O)C(=O)O'),
         ('vanillin',         'COc1cc(C=O)ccc1O'),
         ('catechol',         'Oc1ccccc1O'),
         ('methanol',         'CO'),
         ('ethanol',          'CCO'),
         ('glycolic_acid',    'OCC(=O)O'),
         ('malic_acid',       'OC(=O)CC(O)C(=O)O'),
         ('tartaric_acid',    'OC(=O)C(O)C(O)C(=O)O'),
         ('malonic_acid',     'OC(=O)CC(=O)O'),
         ('ethylene_glycol',  'OCCO'),
         ('benzoic_acid',     'OC(=O)c1ccccc1')],
        ['NH42SO4', 'NH4NO3'],
        [
            # 20 organics each at ~0.01-0.02 mass frac plus 2 salts
            dict(zip(
                ['pinic_acid','pinonic_acid','succinic_acid','glutaric_acid',
                 'adipic_acid','glycerol','glucose','levoglucosan',
                 'acetic_acid','lactic_acid','vanillin','catechol',
                 'methanol','ethanol','glycolic_acid','malic_acid',
                 'tartaric_acid','malonic_acid','ethylene_glycol',
                 'benzoic_acid','NH42SO4','NH4NO3'],
                [0.01]*20 + [0.05, 0.05],
            )),
            dict(zip(
                ['pinic_acid','pinonic_acid','succinic_acid','glutaric_acid',
                 'adipic_acid','glycerol','glucose','levoglucosan',
                 'acetic_acid','lactic_acid','vanillin','catechol',
                 'methanol','ethanol','glycolic_acid','malic_acid',
                 'tartaric_acid','malonic_acid','ethylene_glycol',
                 'benzoic_acid','NH42SO4','NH4NO3'],
                [0.02]*20 + [0.03, 0.03],
            )),
        ],
    ),
]


def decompose(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f'invalid SMILES {smi!r}')
    return {k: v for k, v in groups.aiomfac(mol).items() if v != 0}


def render_input_file(case_id, neutrals, salts_records, composition_table, T_K):
    """neutrals: list of (name, subgroup_dict).
    salts_records: list of (name, cation_id, n_c, anion_id, n_a)
    composition_table: list of dicts {component_name: mass_fraction}"""
    components = [(n, sg) for n, sg in neutrals]
    for name, cid, n_c, aid, n_a in salts_records:
        components.append((name, {cid: int(n_c), aid: int(n_a)}))

    lines = ['Input file for AIOMFAC-web model', '', 'mixture components:', '----']
    for i, (cname, subgroups) in enumerate(components, start=1):
        lines.append(f'component no.:\t{i:02d}')
        lines.append(f"component name:\t'{cname}'")
        for sg_id in sorted(subgroups):
            lines.append(f'subgroup no., qty:\t{sg_id:03d},\t{int(subgroups[sg_id]):02d}')
        lines.append('----')
    lines.append('++++')
    lines.append('mixture composition and temperature:')
    lines.append('mass fraction?\t1')
    lines.append('mole fraction?\t0')
    lines.append('----')
    n_comp = len(components)
    header = ['point', 'T_K'] + [f'cp{i:02d}' for i in range(2, n_comp + 1)]
    lines.append(',\t'.join(header))
    component_names = [c[0] for c in components]
    for k, row in enumerate(composition_table, start=1):
        # row is {component_name -> mass fraction}; water (component 1) is implicit
        vals = [str(k), f'{T_K:.2f}']
        for cname in component_names[1:]:
            vals.append(f'{row.get(cname, 0.0):.5f}')
        lines.append('\t'.join(vals))
    lines.append('====')
    lines.append('')
    return '\n'.join(lines)


def parse_fortran_output(path):
    """Parse the Fortran output. Returns dict with keys:
       n_neutral, x_neutral, gamma_neutral (mole-frac basis a_coeff_x),
       n_cation, n_anion, m_cation, gamma_cation_m, m_anion, gamma_anion_m,
       cation_names, anion_names.
    """
    text = path.read_text()
    # Split into per-species blocks
    blocks = re.split(r"Mixture's component # :|Mixture's species, #  :", text)[1:]
    species = []
    for block in blocks:
        m_name = re.search(r"name\s*:\s*'([^']+)'", block)
        if not m_name:
            continue
        name = m_name.group(1)
        # Determine kind: "Component" -> neutral, "Cation" or "Anion" -> ion
        is_cation = "Cation's" in block
        is_anion = "Anion's" in block
        kind = 'cation' if is_cation else ('anion' if is_anion else 'neutral')
        rows = []
        for line in block.splitlines():
            parts = line.split()
            if len(parts) < 8:
                continue
            try:
                int(parts[0])
            except ValueError:
                continue
            # Columns: no., T, RH, w, x, m, a_coeff_(x|m), a_(x|m), flag
            x = float(parts[4])
            m = float(parts[5])
            gamma = float(parts[6])
            rows.append((x, m, gamma))
        species.append((kind, name, rows))
    return species


def main():
    if not FORTRAN_EXE.exists():
        raise FileNotFoundError(f'Fortran exe not found at {FORTRAN_EXE}')

    fixtures = []
    skipped = []

    with tempfile.TemporaryDirectory(prefix='aiomfac_elec_') as tmp:
        work = Path(tmp)
        (work / 'Inputfiles').mkdir()
        (work / 'Outputfiles').mkdir()

        for ci, (case_id, neutrals, salt_names, comp_table) in enumerate(CASES):
            # Decompose neutrals
            neutral_records = []
            for name, smi in neutrals:
                try:
                    sg = decompose(smi)
                except Exception as e:
                    skipped.append((case_id, str(e)))
                    break
                if any(k in EXTENSION_IDS or k > TOPSUBNO for k in sg):
                    skipped.append((case_id, f'{name} has Fortran-incompatible subgroups {sg}'))
                    break
                neutral_records.append({'name': name, 'smiles': smi, 'subgroups': sg})
            else:
                # Look up salts
                salt_records = []
                for sname in salt_names:
                    matches = [sid for sid, s in data.AIOMFAC_SALT_NAME.items() if s == sname]
                    if not matches:
                        skipped.append((case_id, f'salt {sname!r} not in catalog'))
                        break
                    sid = matches[0]
                    salt_records.append((
                        sname,
                        int(data.AIOMFAC_SALT_CATION_GROUP[sid]),
                        float(data.AIOMFAC_SALT_CATION_STOICH[sid]),
                        int(data.AIOMFAC_SALT_ANION_GROUP[sid]),
                        float(data.AIOMFAC_SALT_ANION_STOICH[sid]),
                    ))
                else:
                    # Build composition rows resolving water mass-fraction
                    rows = []
                    for row in comp_table:
                        total_non_water = sum(row.values())
                        full = {'water': 1.0 - total_non_water, **row}
                        rows.append(full)

                    # Write Fortran input
                    short_id = f'e{ci+1:03d}'
                    input_path = work / 'Inputfiles' / f'input_{short_id}.txt'
                    input_path.write_text(render_input_file(
                        case_id,
                        neutrals=[(r['name'], r['subgroups']) for r in neutral_records],
                        salts_records=salt_records,
                        composition_table=rows,
                        T_K=TEMPERATURE_K,
                    ))
                    # Run Fortran
                    proc = subprocess.run(
                        [str(FORTRAN_EXE), f'./Inputfiles/{input_path.name}'],
                        cwd=work, capture_output=True, text=True, timeout=120,
                    )
                    if proc.returncode != 0:
                        skipped.append((case_id, f'Fortran exit {proc.returncode}: {proc.stderr[:200]}'))
                        continue
                    out_path = work / 'Outputfiles' / f'AIOMFAC_output_{short_id[-4:]}.txt'
                    if not out_path.exists():
                        skipped.append((case_id, f'Fortran output missing at {out_path}'))
                        continue
                    errlog = work / 'Outputfiles' / f'Errorlog_{short_id[-4:]}.txt'
                    if errlog.exists() and 'AIOMFAC ERROR' in errlog.read_text():
                        skipped.append((case_id, f'AIOMFAC error: see {errlog}'))
                        continue
                    species = parse_fortran_output(out_path)
                    fixtures.append({
                        'case_id': case_id,
                        'T_K': TEMPERATURE_K,
                        'neutrals': neutral_records,
                        'salts': salt_names,
                        'mass_fractions_rows': rows,
                        'species': [{'kind': k, 'name': n, 'x_m_gamma': r}
                                    for k, n, r in species],
                    })

    # JSON requires string keys for subgroup dicts
    for fx in fixtures:
        for r in fx['neutrals']:
            r['subgroups'] = {str(k): v for k, v in r['subgroups'].items()}

    FIXTURE_FILE.write_text(json.dumps(fixtures, indent=2) + '\n')
    print(f'Wrote {len(fixtures)} fixtures to {FIXTURE_FILE}')
    if skipped:
        print(f'Skipped {len(skipped)} cases:')
        for cid, reason in skipped:
            print(f'  {cid}: {reason}')


if __name__ == '__main__':
    main()
