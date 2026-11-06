# Third-party material

AeroSysProp is distributed under the GNU General Public License v3 (see
[LICENSE](LICENSE)). It descends from UManSysProp (Topping et al., 2016,
https://doi.org/10.5194/gmd-9-899-2016) and contains material from the
sources below. Each method card in `aerosysprop/methods/cards/` records the
references, parameter files and SHA-256 fingerprints for that method.

## SMARTS pattern sets

| File | Origin | Licence |
|---|---|---|
| `data/FP4.smarts` | Functional group patterns by C. Laggner (Inte:Ligand, 2005), distributed with Open Babel | GNU LGPL. See `data/FP4.smarts.readme` |
| `data/MACCS.smarts` | MACCS key patterns by G. Landrum (RDKit), distributed with Open Babel | BSD. See `data/MACCS.smarts.readme` |

## Group contribution parameters

The coefficients in `aerosysprop/data/` are transcribed from the published
methods listed below, through UManSysProp unless stated otherwise.

| Method | Source of parameters |
|---|---|
| AIOMFAC (SR, MR, LR) | Zuend et al. (2008), https://doi.org/10.5194/acp-8-4559-2008; Zuend et al. (2011), https://doi.org/10.5194/acp-11-9155-2011. Validated against the public Fortran AIOMFAC code (https://github.com/andizuend/AIOMFAC). |
| Joback and Reid | Joback and Reid (1987), https://doi.org/10.1080/00986448708960487 |
| Stein and Brown | Stein and Brown (1994), https://doi.org/10.1021/ci00019a016 |
| Nannoolal boiling point | Nannoolal et al. (2004), https://doi.org/10.1016/j.fluid.2004.09.001 |
| Nannoolal critical properties | Nannoolal et al. (2007), https://doi.org/10.1016/j.fluid.2006.11.014 |
| Nannoolal vapour pressure | Nannoolal et al. (2008), https://doi.org/10.1016/j.fluid.2008.04.020 |
| Myrdal and Yalkowsky | Myrdal and Yalkowsky (1997), https://doi.org/10.1021/ie950242l |
| EVAPORATION | Compernolle et al. (2011), https://doi.org/10.5194/acp-11-9431-2011 |
| SIMPOL.1 | Pankow and Asher (2008), https://doi.org/10.5194/acp-8-2773-2008 |
| Girolami | Girolami (1994), https://doi.org/10.1021/ed071p962 |
| Le Bas, Schroeder, Tyn and Calus | As described in Poling, Prausnitz and O'Connell, The Properties of Gases and Liquids, 5th edn (2001) |

Where AeroSysProp departs from a published method, the departure is recorded
in that method's card under `implementation_notes` and `known_limitations`.

## Optional external tools

OPERA (https://github.com/NIEHS/OPERA) is an optional, separately installed
tool with its own licence. AeroSysProp does not bundle or redistribute it.
See `opera/EXTERNAL_TOOLS.md`.
