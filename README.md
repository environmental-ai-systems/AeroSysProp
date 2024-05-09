# AeroSysProp
A facility to parse chemical structures and predict molecular properties for use in aerosol research.

---

### Badges
![Build Status - just an example for now](https://img.shields.io/badge/build-passing-brightgreen)
![Test Coverage - just an example for now](https://img.shields.io/badge/coverage-90%25-brightgreen)
![License - just an example for now](https://img.shields.io/badge/license-MIT-green)
![Version - just an example for now](https://img.shields.io/badge/version-1.0.0-blue)

## Table of Contents
- [Background](#background)
- [Install](#install)
- [Usage](#usage)
- [Features](#features)
- [Configuration](#configuration)
- [API Documentation](#api-documentation)
- [Contributing](#contributing)
- [Credits and development team](#credits)
- [License](#license)
- [Contact Information](#contact-information)

## Background
### Purpose and Scope
AeroSysProp builds on a previous package, UManSysProp, with the purpose of predicting molecular and mixture properties of relevance to atmospheric aerosol research. These predictions are based on user provided chemical identifiers, such as SMILES strings. These can be individual identifiers or provided as a text file or held within one of multiple database formats. At the heart of AeroSysProp is the reliance on the OpenBabel suite of tools. Using the Python interface to OpenBabel, AeroSysProp connects the many features OpenBabel has to offer with optional suport from other packages in the chemoinformatics ecosystem (including RDKit). Choice of predictive technique, and modifications to any such technique, is determined via the project governance structure. As a minimum, it is our aim to provide a code replica of published methods. Additional modifications to any methods in an attempt to fine-tune, or adapt them, to specific compositions are presented with clear documentation that describe seperating features from the original technique.  

## Setup instructions
### Dependencies
AeroSysProp relies on the following core mdules
- OpenBabel:
- 

### Using the Python code

```bash
# Example command line installation steps
git clone https://yourproject/repository.git
cd repository
npm install
```
### Using the Docker image


## Usage
AeroSysProp is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

## Citation
TBC

## Features

## Configuration

## Documentation

## Contributing

## Credits and development team

## License

## Contact Information

