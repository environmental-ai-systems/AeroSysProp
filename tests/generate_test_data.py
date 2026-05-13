##########################################################################################
# This file generates data to be used with tests. New test data should only be generated when 
# each new method has been documented and agreed on by the management board. This might include
# a review of a parsing scheme or correction of method. Please refer to the guidelines on version
# release governance.
#
# Any new tests that might be use case specific should follow the same pattern provided here
# Each use case that relies on a specific library of SMILES, in addition to the default, should
# be specified in the 
#
# This file will find all functions within the 'methods' folder and generate test data for each one. 
# There is the option to 
#
# Current properties calculated:
# - Pure component vapour pressure
# - Boiling point
# - Subcooled liquid density
# - Critical properties
#
# Last modification 17/05/24

import glob
from rdkit import Chem
import collections
import sys
import os
import json
import dill

parent_directory = os.path.abspath('../..')
sys.path.append(parent_directory)

from aerosysprop.methods import boiling_points
from aerosysprop.methods import vapour_pressures
from aerosysprop.methods import critical_properties
from aerosysprop.methods import liquid_densities

# ------------ This file is seperated into sections --------------
# A) - Individual property (molecule) predictions
# B) - Mixture predictions
#
# Each test dataset is held within a nested dictionary which is then 
# saved to a JSON file for later retrieval. 
# To do this we will need specify each function to be cycle through
# this cannot be done automatically since there is a chance some functions in each 
# module relate to generic capability. 

# 1) Generate dictionary of SMILES objects to pass to functions to generate results

##########################################################################################
# A) -  Individual property prediction
# We can select a number of SMILES string files. Some are sstored as XML files, others plain text
# The following loop works through the files stored in the SMILES folder 


print(‘Mapping species names to SMILES and RDKit mol objects to generate test data’)

smiles_folder = ‘/SMILES’

SMILES_dict=collections.defaultdict(lambda: collections.defaultdict)
Pybel_object_dict=collections.defaultdict(lambda: collections.defaultdict)

print(‘Mapping species names to SMILES and RDKit mol objects’)

for filename in os.listdir(smiles_folder):

    if filename.endswith(‘.txt’):
        text=open(filename[:],’r’)
        for line in text:
            input = line.split()
            try:
                species_SMILES=input[1]
                SMILES_dict[filename][species_SMILES]=species_SMILES
                mol = Chem.MolFromSmiles(input[1])
                if mol is None:
                    raise ValueError(‘Could not parse SMILES: %s’ % input[1])
                Pybel_object_dict[filename][species_SMILES]=mol
            except:
                print(‘Problem generating mol object from SMILES string’)
    elif filename.endswith(‘.xml’):
        dom = minidom.parse(SMILES_filename)
        for node in dom.getElementsByTagName(‘species’):
            species_name=node.attributes[‘species_name’].value
            if species_name not in [‘OH’,’HO2’,’O3’,’hv’]:
                if species_name in species_dict.values():
                    try:
                        species_SMILES=re.search(r’<smiles>(.*?)</smiles>’,node.getElementsByTagName(‘smiles’).item(0).toxml()).group(1)
                        SMILES_dict[filename][species_name]=species_SMILES
                        mol = Chem.MolFromSmiles(species_SMILES)
                        if mol is None:
                            raise ValueError(‘Could not parse SMILES: %s’ % species_SMILES)
                        Pybel_object_dict[filename][species_name]=mol
                    except:
                        print(‘No SMILES entry for species ‘,species_name)

with open(‘SMILES_master.json’, ‘w’) as fp:
    json.dump(SMILES_dict, fp)

# Now generated the test data from the collected SMILES strings.

step=0
filenames=[]
Compound_reference=[]
smiles_array=[]
property_array=[]
Pybel_object_dict=dict()

for filename in os.listdir(smiles_folder):
    
    # 2a) Boiling points [(K)]
    boiling_point_dict=collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))

    # Now leed to loop through 'Pybel_object_dict[filename][species_name]=Pybel_object'

    for species_name, Pybel_object in Pybel_object_dict[filename].items():
        
        boiling_point_dict[filename][species_name]['joback_and_reid']=boiling_points.joback_and_reid(Pybel_object)
        boiling_point_dict[filename][species_name]['stein_and_brown']=boiling_points.stein_and_brown(Pybel_object)
        boiling_point_dict[filename][species_name]['nannoolal']=boiling_points.nannoolal(Pybel_object)
    
