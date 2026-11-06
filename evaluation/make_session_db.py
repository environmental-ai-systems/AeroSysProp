"""Write the session database the evaluation's read_session_database questions use.

    python -m evaluation.make_session_db

Writes evaluation/data/session.db the way the web app saves a session: the
results of batch.run_all for a few SOA-relevant compounds, at 298.15 K and at
273.15 K. Its values come from the installed methods, the same as the keys, so
re-run this (and then build_questions) after any approved parameter change.

The two temperatures matter: pinic acid is less volatile than 2-methyltetrol
by SIMPOL at 298.15 K and the other way round at 273.15 K, so an answer that
reads the wrong temperature's rows gives the wrong compound.
"""
from __future__ import annotations

import sys
from pathlib import Path

from rdkit import Chem

from aerosysprop import batch, db

HERE = Path(__file__).resolve().parent
OUT = HERE / 'data' / 'session.db'

COMPOUNDS = [
    ('malonic acid', 'OC(=O)CC(=O)O'),
    ('succinic acid', 'OC(=O)CCC(=O)O'),
    ('glutaric acid', 'OC(=O)CCCC(=O)O'),
    ('adipic acid', 'OC(=O)CCCCC(=O)O'),
    ('pinonic acid', 'CC(=O)C1CC(C1(C)C)CC(=O)O'),
    ('pinic acid', 'CC1(C)C(CC(=O)O)CC1C(=O)O'),
    ('2-methyltetrol', 'CC(O)(CO)C(O)CO'),
    ('levoglucosan', 'OC1C(O)C2COC(O2)C(O)C1O'),
]
METHODS = ['vp_simpol', 'vp_evaporation', 'vp_nannoolal']
TEMPERATURES = [298.15, 273.15]


def main() -> int:
    OUT.unlink(missing_ok=True)
    mols = [(name, Chem.MolFromSmiles(smi)) for name, smi in COMPOUNDS]
    for T in TEMPERATURES:
        results = batch.run_all(mols, methods=METHODS, temperature=T, use_threads=True)
        db.save_results(str(OUT), results, T)
    print(f'{len(COMPOUNDS)} compounds at {TEMPERATURES} K -> {OUT}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
