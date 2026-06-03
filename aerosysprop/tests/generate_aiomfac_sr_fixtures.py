"""Generate Fortran-AIOMFAC reference outputs for the SR regression tests.

For each entry in CURATED_BINARIES and MULTI_COMPONENT, this script:
  1. Decomposes each SMILES via ``groups.aiomfac()``.
  2. Filters out compounds whose decomposition contains subgroup IDs the
     Fortran reference can't handle (>265, the 281-285 Mark-Barley extensions,
     fractional / negative counts).
  3. Writes an AIOMFAC input file at a chosen sweep of compositions.
  4. Runs the local Fortran build (../AIOMFAC/FortranCode/AIOMFAC-web.out).
  5. Parses the ``a_coeff_x`` column from the output text file.
  6. Aggregates everything into ``aiomfac_sr_fixtures.json`` alongside this
     script.

The Fortran executable is expected at
``/Users/davidtopping/Documents/GitHub/AIOMFAC/FortranCode/AIOMFAC-web.out``.
Set ``SDKROOT=$(xcrun --show-sdk-path)`` first if invoking on macOS so the
linker (and thereby any rebuild) finds libSystem.

Re-run this script whenever the curated list, the AIOMFAC parameters, or the
``groups.aiomfac()`` matcher changes.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from aerosysprop import groups


HERE = Path(__file__).parent
FIXTURE_FILE = HERE / 'aiomfac_sr_fixtures.json'
SKIPPED_LOG = HERE / 'aiomfac_sr_skipped.txt'

FORTRAN_DIR = Path('/Users/davidtopping/Documents/GitHub/AIOMFAC/FortranCode')
FORTRAN_EXE = FORTRAN_DIR / 'AIOMFAC-web.out'

# Fortran AIOMFAC limit on subgroup IDs (ModSystemProp.f90:54).
TOPSUBNO = 265

# Subgroups produced by aerosysprop's "Mark Barley 2012 extensions" that the
# Fortran reference will silently drop (out of range).
EXTENSION_IDS = {281, 282, 283, 284, 285}

# Water as an AIOMFAC subgroup dict (subgroup 16 is the special H2O group).
WATER_SUBGROUPS = {16: 1}

# Mole-fraction sweep applied to every binary (x_organic).
BINARY_X_SWEEP = (0.01, 0.05, 0.10, 0.30, 0.50, 0.70, 0.90, 0.95, 0.99)

# Curated list of binary water + organic test compounds.
CURATED_BINARIES = [
    # ---- alcohols
    ('methanol',         'CO'),
    ('ethanol',          'CCO'),
    ('1-propanol',       'CCCO'),
    ('2-propanol',       'CC(O)C'),
    ('1-butanol',        'CCCCO'),
    ('tert-butanol',     'CC(C)(C)O'),
    # ---- diols / polyols
    ('ethylene_glycol',  'OCCO'),
    ('1,4-butanediol',   'OCCCCO'),
    ('glycerol',         'OCC(O)CO'),
    ('erythritol',       'OCC(O)C(O)CO'),
    ('xylitol',          'OCC(O)C(O)C(O)CO'),
    # ---- carboxylic acids
    ('formic_acid',      'OC=O'),
    ('acetic_acid',      'CC(=O)O'),
    ('propionic_acid',   'CCC(=O)O'),
    ('oxalic_acid',      'OC(=O)C(=O)O'),
    ('malonic_acid',     'OC(=O)CC(=O)O'),
    ('succinic_acid',    'OC(=O)CCC(=O)O'),
    ('glutaric_acid',    'OC(=O)CCCC(=O)O'),
    ('adipic_acid',      'OC(=O)CCCCC(=O)O'),
    ('benzoic_acid',     'OC(=O)c1ccccc1'),
    # ---- hydroxy / multifunctional acids
    ('lactic_acid',      'CC(O)C(=O)O'),
    ('glycolic_acid',    'OCC(=O)O'),
    ('malic_acid',       'OC(=O)CC(O)C(=O)O'),
    ('tartaric_acid',    'OC(=O)C(O)C(O)C(=O)O'),
    ('citric_acid',      'OC(=O)CC(O)(C(=O)O)CC(=O)O'),
    ('pyruvic_acid',     'CC(=O)C(=O)O'),
    # ---- atmospheric (alpha-pinene oxidation products)
    ('pinic_acid',       'CC1(C)C(C(=O)O)CCC1CC(=O)O'),
    ('pinonic_acid',     'CC(=O)C1CCC(C1(C)C)CC(=O)O'),
    # ---- aldehydes
    ('formaldehyde',     'C=O'),
    ('acetaldehyde',     'CC=O'),
    ('propanal',         'CCC=O'),
    # ---- ketones
    ('MEK',              'CCC(=O)C'),
    # ---- esters
    ('methyl_acetate',   'COC(=O)C'),
    ('ethyl_acetate',    'CCOC(=O)C'),
    # ---- ethers
    ('diethyl_ether',    'CCOCC'),
    ('THF',              'C1CCOC1'),
    # ---- PEG-like
    ('diethylene_glycol','OCCOCCO'),
    # ---- aromatics
    ('benzene',          'c1ccccc1'),
    ('toluene',          'Cc1ccccc1'),
    ('xylene',           'Cc1ccccc1C'),
    ('phenol',           'Oc1ccccc1'),
    ('m-cresol',         'Cc1cccc(O)c1'),
    ('catechol',         'Oc1ccccc1O'),
    ('resorcinol',       'Oc1cccc(O)c1'),
    ('benzaldehyde',     'O=Cc1ccccc1'),
    ('vanillin',         'COc1cc(C=O)ccc1O'),
    # ---- sugars
    ('levoglucosan',     'OC[C@H]1OC2OC1C(O)C2O'),
    ('glucose',          'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'),
    ('sucrose',          'OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O'),
    # ---- alkenes / alkynes
    ('1-butene',         'CCC=C'),
    ('cis-2-butene',     'C/C=C/C'),
    ('isoprene',         'CC(=C)C=C'),
    # ---- amines (matcher's heteroatom coverage is partial; mismatches will
    # be diagnosed in the generator log)
    ('methylamine',      'CN'),
    ('ethylamine',       'CCN'),
    ('dimethylamine',    'CNC'),
    ('trimethylamine',   'CN(C)C'),
    ('ethanolamine',     'OCCN'),
    ('aniline',          'c1ccc(N)cc1'),
    # ---- amides
    ('acetamide',        'CC(=O)N'),
    ('DMF',              'CN(C)C=O'),
    # ---- nitriles
    ('acetonitrile',     'CC#N'),
    # ---- halides
    ('chloroform',       'ClC(Cl)Cl'),
    ('CCl4',             'ClC(Cl)(Cl)Cl'),
    ('chlorobenzene',    'Clc1ccccc1'),
    # ---- hydrocarbons (extended)
    ('n-hexane',         'CCCCCC'),
    ('octane',           'CCCCCCCC'),
    ('cyclohexane',      'C1CCCCC1'),
    ('methylcyclohexane','CC1CCCCC1'),
    # ---- atmospheric oxidation products from tests/SMILES/short_list.txt
    # (those that pre-filter as Fortran-compatible)
    ('SL_02',            'CC1(C2CCC(O)C1C2)C'),
    ('SL_03',            'C(C(O)OON(O)O)C1C(C(C(O)C)C1)(C)C'),
    ('SL_05',            'C(O)(OON(O)O)CO'),
    ('SL_06',            'C1(C(C(CCO)C1)(C)C)C(O)C'),
    ('SL_09',            'CC1(C2CCC(CO)(C1C2)ON(O)O)C'),
    ('SL_11',            'OCCC(O)OON(O)O'),
    ('SL_14',            'C(O)CO'),
    ('SL_15',            'CC(C(O)OON(O)O)C'),
    ('SL_16',            'CC(O)C(O)CC(O)OON(O)O'),
    ('SL_17',            'C12C(C(CC(C1(C)O)ON(O)O)C2)(C)C'),
    ('SL_18',            'C(C(O)O)C1C(C(C(C)O)C1)(C)C'),
    ('SL_19',            'CC(C)(C(O)OON(O)O)O'),
]

# Multi-component mixtures: list of (mix_name, [(component_name, smiles), ...],
# [composition_dict, ...]) where each composition_dict maps component_name -> x.
# The first composition entry is for water; component-name "water" is implicit.
MULTI_COMPONENT = [
    (
        'water_ethanol_ethyleneglycol',
        [('ethanol', 'CCO'), ('ethylene_glycol', 'OCCO')],
        [
            {'ethanol': 0.10, 'ethylene_glycol': 0.10},
            {'ethanol': 0.30, 'ethylene_glycol': 0.30},
            {'ethanol': 0.50, 'ethylene_glycol': 0.20},
            {'ethanol': 0.10, 'ethylene_glycol': 0.80},
            {'ethanol': 0.80, 'ethylene_glycol': 0.10},
        ],
    ),
    (
        'water_aceticacid_ethylacetate',
        [('acetic_acid', 'CC(=O)O'), ('ethyl_acetate', 'CCOC(=O)C')],
        [
            {'acetic_acid': 0.10, 'ethyl_acetate': 0.10},
            {'acetic_acid': 0.30, 'ethyl_acetate': 0.30},
            {'acetic_acid': 0.50, 'ethyl_acetate': 0.20},
            {'acetic_acid': 0.20, 'ethyl_acetate': 0.50},
        ],
    ),
    (
        'water_phenol_toluene',
        [('phenol', 'Oc1ccccc1'), ('toluene', 'Cc1ccccc1')],
        [
            {'phenol': 0.10, 'toluene': 0.10},
            {'phenol': 0.30, 'toluene': 0.10},
            {'phenol': 0.10, 'toluene': 0.30},
            {'phenol': 0.40, 'toluene': 0.40},
        ],
    ),
    (
        'water_glucose_aceticacid',
        [('glucose', 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'),
         ('acetic_acid', 'CC(=O)O')],
        [
            {'glucose': 0.05, 'acetic_acid': 0.10},
            {'glucose': 0.10, 'acetic_acid': 0.30},
            {'glucose': 0.20, 'acetic_acid': 0.10},
        ],
    ),
    (
        'water_ethanol_aceticacid_glucose',
        [('ethanol', 'CCO'), ('acetic_acid', 'CC(=O)O'),
         ('glucose', 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O')],
        [
            {'ethanol': 0.20, 'acetic_acid': 0.20, 'glucose': 0.10},
            {'ethanol': 0.40, 'acetic_acid': 0.20, 'glucose': 0.05},
            {'ethanol': 0.10, 'acetic_acid': 0.10, 'glucose': 0.20},
        ],
    ),
    # mixed alcohol + diol + ester
    (
        'water_ethanol_ethyleneglycol_ethylacetate',
        [('ethanol', 'CCO'), ('ethylene_glycol', 'OCCO'),
         ('ethyl_acetate', 'CCOC(=O)C')],
        [
            {'ethanol': 0.20, 'ethylene_glycol': 0.20, 'ethyl_acetate': 0.10},
            {'ethanol': 0.10, 'ethylene_glycol': 0.30, 'ethyl_acetate': 0.30},
            {'ethanol': 0.30, 'ethylene_glycol': 0.10, 'ethyl_acetate': 0.40},
        ],
    ),
    # alcohol + ether + aromatic
    (
        'water_ethanol_THF_toluene',
        [('ethanol', 'CCO'), ('THF', 'C1CCOC1'), ('toluene', 'Cc1ccccc1')],
        [
            {'ethanol': 0.20, 'THF': 0.20, 'toluene': 0.10},
            {'ethanol': 0.10, 'THF': 0.30, 'toluene': 0.10},
            {'ethanol': 0.40, 'THF': 0.20, 'toluene': 0.20},
        ],
    ),
    # dicarboxylic acids + alcohol (atmospheric aging proxy)
    (
        'water_succinic_glutaric_ethanol',
        [('succinic_acid', 'OC(=O)CCC(=O)O'),
         ('glutaric_acid', 'OC(=O)CCCC(=O)O'),
         ('ethanol',       'CCO')],
        [
            {'succinic_acid': 0.05, 'glutaric_acid': 0.05, 'ethanol': 0.10},
            {'succinic_acid': 0.10, 'glutaric_acid': 0.10, 'ethanol': 0.20},
            {'succinic_acid': 0.15, 'glutaric_acid': 0.05, 'ethanol': 0.30},
        ],
    ),
    # multi-functional atmospheric oxidation product mix
    (
        'water_pinicacid_pinonicacid_levoglucosan',
        [('pinic_acid',   'CC1(C)C(C(=O)O)CCC1CC(=O)O'),
         ('pinonic_acid', 'CC(=O)C1CCC(C1(C)C)CC(=O)O'),
         ('levoglucosan', 'OC[C@H]1OC2OC1C(O)C2O')],
        [
            {'pinic_acid': 0.05, 'pinonic_acid': 0.05, 'levoglucosan': 0.05},
            {'pinic_acid': 0.10, 'pinonic_acid': 0.05, 'levoglucosan': 0.10},
            {'pinic_acid': 0.05, 'pinonic_acid': 0.10, 'levoglucosan': 0.15},
        ],
    ),
]

TEMPERATURE_K = 298.15

# MCM mechanism file: we dynamically pull every closed-shell species with a
# unique AIOMFAC-subgroup signature as a water-binary fixture. Radicals (any
# bracketed atom token like [O], [CH2], etc.) are excluded since they aren't
# stable species and Fortran AIOMFAC isn't applicable to them.
MCM_XML = HERE / 'SMILES' / 'MCM331.xml'
INCLUDE_MCM = True
# Cap on the number of MCM unique-signature binaries to add (set to a small
# number for fast iteration; None for everything).
MCM_LIMIT = None


def decompose(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'invalid SMILES {smiles!r}')
    return {k: v for k, v in groups.aiomfac(mol).items() if v != 0}


def check_fortran_compatible(name, subgroups):
    """Return None if ok, else a reason string for skipping."""
    out_of_range = [k for k in subgroups if k > TOPSUBNO or k < 1]
    if out_of_range:
        return f'subgroup IDs out of Fortran range: {sorted(out_of_range)}'
    extensions = [k for k in subgroups if k in EXTENSION_IDS]
    if extensions:
        return f'aerosysprop extension subgroups 281-285 (no Fortran counterpart): {sorted(extensions)}'
    fractional = [k for k, v in subgroups.items() if abs(v - round(v)) > 1e-9]
    if fractional:
        return f'fractional subgroup counts: {[(k, subgroups[k]) for k in fractional]}'
    negative = [k for k, v in subgroups.items() if v < 0]
    if negative:
        return f'negative subgroup counts: {[(k, subgroups[k]) for k in negative]}'
    return None


def render_input_file(case_id, components, composition_table, T_K=TEMPERATURE_K):
    """Build an AIOMFAC text input. ``components`` is a list of
    (display_name, {subgroup_id: int_count}) starting with water. The
    composition_table is a list of dicts giving cp02, cp03, ... values
    keyed by 1-based component index (cp02 means component 2's x)."""
    lines = []
    lines.append('Input file for AIOMFAC-web model')
    lines.append('')
    lines.append('mixture components:')
    lines.append('----')
    for i, (cname, subgroups) in enumerate(components, start=1):
        lines.append(f'component no.:\t{i:02d}')
        lines.append(f"component name:\t'{cname}'")
        for sg_id in sorted(subgroups):
            count = int(subgroups[sg_id])
            lines.append(f'subgroup no., qty:\t{sg_id:03d},\t{count:02d}')
        lines.append('----')
    lines.append('++++')
    lines.append('mixture composition and temperature:')
    lines.append('mass fraction?\t0')
    lines.append('mole fraction?\t1')
    lines.append('----')
    # composition header: "point, T_K, cp02, cp03, ..."
    header = ['point', 'T_K'] + [f'cp{i:02d}' for i in range(2, len(components) + 1)]
    lines.append(',\t'.join(header))
    for k, row in enumerate(composition_table, start=1):
        vals = [str(k), f'{T_K:.2f}'] + [
            f'{row[i]:.5f}' for i in range(2, len(components) + 1)
        ]
        lines.append('\t'.join(vals))
    lines.append('====')
    lines.append('')
    return '\n'.join(lines)


def parse_fortran_output(path):
    """Return (n_points x n_components) arrays of x_i and gamma_x_i parsed from
    the Fortran AIOMFAC output text file."""
    text = path.read_text()
    blocks = re.split(r"Mixture's component # :", text)[1:]
    parsed = []
    for block in blocks:
        rows = []
        for line in block.splitlines():
            parts = line.split()
            if len(parts) < 8:
                continue
            try:
                int(parts[0])
            except ValueError:
                continue
            x = float(parts[4])
            gamma = float(parts[6])
            rows.append((x, gamma))
        parsed.append(rows)
    n_pts = len(parsed[0])
    n_comp = len(parsed)
    assert all(len(r) == n_pts for r in parsed), 'component blocks misaligned'
    x_matrix = [[parsed[j][k][0] for j in range(n_comp)] for k in range(n_pts)]
    gamma_matrix = [[parsed[j][k][1] for j in range(n_comp)] for k in range(n_pts)]
    return x_matrix, gamma_matrix


def run_fortran(input_path, work_dir):
    """Run AIOMFAC-web.out on input_path. ``work_dir`` will become the working
    directory; it must contain an ``Inputfiles/`` subdirectory holding the
    input, and an ``Outputfiles/`` subdirectory will receive results."""
    proc = subprocess.run(
        [str(FORTRAN_EXE), f'./Inputfiles/{input_path.name}'],
        cwd=work_dir, capture_output=True, text=True, timeout=120,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f'AIOMFAC executable returned {proc.returncode}\n'
            f'stdout: {proc.stdout}\nstderr: {proc.stderr}'
        )
    # Map back to output file: AIOMFAC uses LAST 4 chars of input basename.
    stem = input_path.stem
    out_path = work_dir / 'Outputfiles' / f'AIOMFAC_output_{stem[-4:]}.txt'
    if not out_path.exists():
        raise FileNotFoundError(f'expected output at {out_path}, not found')
    # Confirm no error from AIOMFAC's own error log
    errlog = work_dir / 'Outputfiles' / f'Errorlog_{stem[-4:]}.txt'
    if errlog.exists():
        text = errlog.read_text()
        if 'AIOMFAC ERROR' in text:
            raise RuntimeError(f'AIOMFAC error in {errlog}: {text}')
    return out_path


def make_binary_compositions(x_organic_values):
    """For a single water+organic binary, build composition_table entries with
    cp02 = x_organic (water is implicit cp01 = 1 - cp02)."""
    return [{2: x} for x in x_organic_values]


def make_multi_compositions(component_names, composition_specs):
    """For multi-component, translate {name: x} dicts to cp02, cp03, ..."""
    table = []
    for spec in composition_specs:
        row = {}
        for j, name in enumerate(component_names, start=2):
            row[j] = spec[name]
        table.append(row)
    return table


def _load_mcm_unique_signatures():
    """Parse MCM XML, return list of (name, smiles, subgroup_dict) with one
    representative per unique subgroup signature. Skips radicals, empty
    decompositions, and Fortran-incompatible subgroup combinations."""
    import re as _re
    if not MCM_XML.exists():
        return []
    text = MCM_XML.read_text()
    pattern = _re.compile(
        r'<species\s+species_number="[^"]+"\s+species_name="([^"]+)"\s*>\s*<smiles>([^<]+)</smiles>',
        _re.DOTALL,
    )
    species = pattern.findall(text)
    sigs = {}
    for name, smi in species:
        # Skip anything with bracketed atom tokens (radicals like [O],
        # carbenes [CH2], etc.). Allow charged ions ([O-], [N+]) as they
        # appear in some nitro/peroxynitrate compounds.
        if _re.search(r'\[[A-Za-z][A-Za-z]?\]', smi):
            continue
        try:
            g = decompose(smi)
        except ValueError:
            continue
        if not g:
            continue
        if check_fortran_compatible(name, g) is not None:
            continue
        sig = frozenset(g.items())
        if sig not in sigs:
            sigs[sig] = (name, smi, g)
    out = list(sigs.values())
    out.sort(key=lambda r: r[0])  # deterministic ordering
    return out


def generate_fixtures():
    if not FORTRAN_EXE.exists():
        raise FileNotFoundError(
            f'Fortran AIOMFAC executable not found at {FORTRAN_EXE}; '
            'build it first (see Fortran build commands).'
        )

    fixtures = []
    skipped = []

    # Use a private working dir so we don't clobber any local outputs the user
    # might have under FORTRAN_DIR/Outputfiles.
    with tempfile.TemporaryDirectory(prefix='aiomfac_fixtures_') as tmp:
        work = Path(tmp)
        (work / 'Inputfiles').mkdir()
        (work / 'Outputfiles').mkdir()

        case_idx = 0

        # ---- Binaries -------------------------------------------------------
        for name, smi in CURATED_BINARIES:
            try:
                sg = decompose(smi)
            except ValueError as e:
                skipped.append((name, smi, str(e)))
                continue
            reason = check_fortran_compatible(name, sg)
            if reason:
                skipped.append((name, smi, reason))
                continue

            case_idx += 1
            short_id = f'b{case_idx:03d}'
            case_id = f'binary__{name}'
            input_name = f'input_{short_id}.txt'
            input_path = work / 'Inputfiles' / input_name
            input_path.write_text(render_input_file(
                case_id,
                components=[('Water', WATER_SUBGROUPS), (name, sg)],
                composition_table=make_binary_compositions(BINARY_X_SWEEP),
            ))
            out_path = run_fortran(input_path, work)
            x_matrix, gamma_matrix = parse_fortran_output(out_path)
            fixtures.append({
                'case_id': case_id,
                'kind': 'binary',
                'T_K': TEMPERATURE_K,
                'components': [
                    {'name': 'water', 'smiles': 'O', 'subgroups': WATER_SUBGROUPS},
                    {'name': name, 'smiles': smi, 'subgroups': sg},
                ],
                'x':     x_matrix,
                'gamma': gamma_matrix,
            })

        # ---- MCM unique-signature binaries ----------------------------------
        if INCLUDE_MCM:
            mcm = _load_mcm_unique_signatures()
            if MCM_LIMIT is not None:
                mcm = mcm[:MCM_LIMIT]
            print(f'MCM331: {len(mcm)} unique-signature water-binary cases')
            for name, smi, sg in mcm:
                case_idx += 1
                short_id = f'x{case_idx:04d}'
                case_id = f'mcm__{name}'
                input_name = f'input_{short_id}.txt'
                input_path = work / 'Inputfiles' / input_name
                input_path.write_text(render_input_file(
                    case_id,
                    components=[('Water', WATER_SUBGROUPS), (name, sg)],
                    composition_table=make_binary_compositions(BINARY_X_SWEEP),
                ))
                try:
                    out_path = run_fortran(input_path, work)
                except RuntimeError as e:
                    skipped.append((f'MCM::{name}', smi, f'Fortran error: {e}'[:200]))
                    continue
                x_matrix, gamma_matrix = parse_fortran_output(out_path)
                fixtures.append({
                    'case_id': case_id,
                    'kind': 'mcm_binary',
                    'T_K': TEMPERATURE_K,
                    'components': [
                        {'name': 'water', 'smiles': 'O', 'subgroups': WATER_SUBGROUPS},
                        {'name': name, 'smiles': smi, 'subgroups': sg},
                    ],
                    'x':     x_matrix,
                    'gamma': gamma_matrix,
                })

        # ---- Multi-component ------------------------------------------------
        for mix_name, comp_specs, composition_specs in MULTI_COMPONENT:
            comp_records = []
            usable = True
            for cname, smi in comp_specs:
                try:
                    sg = decompose(smi)
                except ValueError as e:
                    skipped.append((f'{mix_name}::{cname}', smi, str(e)))
                    usable = False
                    break
                reason = check_fortran_compatible(cname, sg)
                if reason:
                    skipped.append((f'{mix_name}::{cname}', smi, reason))
                    usable = False
                    break
                comp_records.append({'name': cname, 'smiles': smi, 'subgroups': sg})
            if not usable:
                continue

            case_idx += 1
            short_id = f'm{case_idx:03d}'
            case_id = f'multi__{mix_name}'
            input_name = f'input_{short_id}.txt'
            input_path = work / 'Inputfiles' / input_name

            input_components = [('Water', WATER_SUBGROUPS)] + [
                (r['name'], r['subgroups']) for r in comp_records
            ]
            input_path.write_text(render_input_file(
                case_id,
                components=input_components,
                composition_table=make_multi_compositions(
                    [r['name'] for r in comp_records], composition_specs,
                ),
            ))
            out_path = run_fortran(input_path, work)
            x_matrix, gamma_matrix = parse_fortran_output(out_path)
            fixtures.append({
                'case_id': case_id,
                'kind': 'multi',
                'T_K': TEMPERATURE_K,
                'components': [
                    {'name': 'water', 'smiles': 'O', 'subgroups': WATER_SUBGROUPS}
                ] + comp_records,
                'x':     x_matrix,
                'gamma': gamma_matrix,
            })

    # Serialise. JSON requires string keys, so encode subgroup IDs as strings.
    def _stringify_subgroups(rec):
        for c in rec['components']:
            c['subgroups'] = {str(k): v for k, v in c['subgroups'].items()}
        return rec

    fixtures = [_stringify_subgroups(r) for r in fixtures]
    FIXTURE_FILE.write_text(json.dumps(fixtures, indent=2) + '\n')

    SKIPPED_LOG.write_text(
        'AIOMFAC SR fixture generation - skipped compounds\n'
        '=================================================\n\n'
        + ''.join(
            f'  {name:40s}  {smi:55s}  {reason}\n' for name, smi, reason in skipped
        )
    )

    print(f'Wrote {len(fixtures)} fixtures to {FIXTURE_FILE}')
    print(f'Skipped {len(skipped)} compounds (see {SKIPPED_LOG})')
    if skipped:
        print('Skipped:')
        for name, smi, reason in skipped:
            print(f'  {name:40s}  {reason}')


if __name__ == '__main__':
    generate_fixtures()
