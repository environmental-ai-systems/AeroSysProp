"""Generate evaporation_set_fixtures.json: MCM species on which EVAPORATION
parameter sets 1 and 2 give different vapour pressures.

Run from this directory only when a parameter set is deliberately changed
(which needs core team approval; see METHODS_GOVERNANCE.md):

    python generate_evaporation_set_fixtures.py
"""
import json
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path

from rdkit import Chem, RDLogger

from aerosysprop import groups
from aerosysprop.methods import vapour_pressures as vp

HERE = Path(__file__).parent
MCM_XML = HERE / 'SMILES' / 'MCM331.xml'
OUT = HERE / 'evaporation_set_fixtures.json'
T = 298.15
ENOL = Chem.MolFromSmarts('[OX2H][CX3]=[CX3]')
ALPHA_HYDROXY_ACID = Chem.MolFromSmarts('[OX2H][CX4][CX3](=O)[OX2H]')
# Controls: molecules with alcohol or acid groups where both sets must agree.
CONTROLS = ['CCO', 'OCCO', 'CC(=O)O', 'OC(=O)CC(=O)O', 'OCC(O)CO', 'CC(O)CC(=O)O', 'C=CCO']


def mcm_species():
    ns = {'m': 'http://mcm.leeds.ac.uk/MCM/mechanism.xsd'}
    root = ET.parse(MCM_XML).getroot()
    for sp in root.iter('{%s}species' % ns['m']):
        smi = sp.find('m:smiles', ns)
        if smi is not None and smi.text:
            yield sp.get('species_name'), smi.text.strip()


def nonzero(d):
    return {k: v for k, v in d.items() if v}


def category(mol):
    if mol.HasSubstructMatch(ENOL):
        return 'enol'
    if mol.HasSubstructMatch(ALPHA_HYDROXY_ACID):
        return 'alpha_hydroxy_acid'
    return 'other'


def main():
    RDLogger.DisableLog('rdApp.*')
    warnings.simplefilter('ignore')
    differing, checked = [], 0
    for name, smi in mcm_species():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        try:
            v1, v2 = vp.evaporation(mol, T), vp.evaporation2(mol, T)
        except Exception:
            continue
        checked += 1
        if abs(v1 - v2) > 1e-9:
            g1, g2 = nonzero(groups.evaporation(mol)), nonzero(groups.evaporation2(mol))
            differing.append({
                'name': name, 'smiles': smi, 'category': category(mol),
                'set1_log10atm': v1, 'set2_log10atm': v2,
                'group_differences': {k: [g1.get(k, 0), g2.get(k, 0)]
                                      for k in sorted(set(g1) | set(g2))
                                      if g1.get(k, 0) != g2.get(k, 0)}})
    controls = []
    for smi in CONTROLS:
        mol = Chem.MolFromSmiles(smi)
        controls.append({'smiles': smi, 'log10atm': vp.evaporation(mol, T)})
    out = {'temperature_K': T, 'source': 'MCM v3.3.1 (tests/SMILES/MCM331.xml)',
           'species_checked': checked, 'differing': differing, 'controls': controls}
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    cats = {}
    for d in differing:
        cats[d['category']] = cats.get(d['category'], 0) + 1
    print(f'{checked} species checked, {len(differing)} differ: {cats}')


if __name__ == '__main__':
    main()
