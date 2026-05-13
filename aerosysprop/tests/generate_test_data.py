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
# WARNING - running this file will overwrite any existing test data. Please ensure you archive, where
# appropriate, older test outputs.
#
# Last modification 17/05/24

import glob
from rdkit import Chem
import collections
import sys
import os
import json
import pickle
from xml.dom import minidom
import re
import tqdm
import numpy as np

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


print('Mapping species names to SMILES and RDKit mol objects to generate test data')

smiles_folder = os.getcwd()+'/SMILES'

SMILES_dict=collections.defaultdict(lambda: collections.defaultdict())
Pybel_object_dict=collections.defaultdict(lambda: collections.defaultdict())

print('Mapping species names to SMILES and RDKit mol objects')

for filename in os.listdir(smiles_folder):

    if filename.endswith('.txt'):
        text=open(smiles_folder+'/'+filename,'r')
        for line in text:
            input = line.split()
            try:
                species_SMILES=input[1]
                SMILES_dict[filename][species_SMILES]=species_SMILES
                mol = Chem.MolFromSmiles(input[1])
                if mol is None:
                    raise ValueError('Could not parse SMILES: %s' % input[1])
                Pybel_object_dict[filename][species_SMILES]=mol
            except:
                print('Problem generating mol object from SMILES string')
    elif filename.endswith('.xml'):
        dom = minidom.parse(smiles_folder+'/'+filename)
        for node in dom.getElementsByTagName('species'):
            species_name=node.attributes['species_name'].value
            if species_name not in ['OH','HO2','O3','hv']:
                try:
                    species_SMILES=re.search(r'<smiles>(.*?)</smiles>',node.getElementsByTagName('smiles').item(0).toxml()).group(1)
                    SMILES_dict[filename][species_name]=species_SMILES
                    mol = Chem.MolFromSmiles(species_SMILES)
                    if mol is None:
                        raise ValueError('Could not parse SMILES: %s' % species_SMILES)
                    Pybel_object_dict[filename][species_name]=mol
                except:
                    print('No SMILES entry for species ',species_name)

# Save the SMILES dictionarites to a JSON file for later loading                        
with open("SMILES_master.json", "w") as fp:
          json.dump(SMILES_dict, fp)
# It seems too difficult to pickle the Pybel object so going to generate it each time


# Now generated the test data from the collected SMILES strings.

step=0
filenames=[]
Compound_reference=[]
smiles_array=[]
property_array=[]

print('Predicting properties and generating test data')

boiling_point_dict=collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict())))
vapour_pressure_dict_BP=collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))))
vapour_pressure_dict=collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict())))
temperature=298.15
density_CP=collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))))
density=collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict())))


for filename in os.listdir(smiles_folder):
    
    # Now leed to loop through 'Pybel_object_dict[filename][species_name]=Pybel_object'
    print('Boiling points, file:',filename)

    for species_name, Pybel_object in tqdm.tqdm(Pybel_object_dict[filename].items()):

        boiling_point_dict[filename][species_name]['joback_and_reid']=boiling_points.joback_and_reid(Pybel_object)
        boiling_point_dict[filename][species_name]['stein_and_brown']=boiling_points.stein_and_brown(Pybel_object)
        boiling_point_dict[filename][species_name]['nannoolal']=boiling_points.nannoolal(Pybel_object)
        
    # Now leed to loop through 'Pybel_object_dict[filename][species_name]=Pybel_object'
    print('Vapour pressure (psat), file:',filename)


    for species_name, Pybel_object in tqdm.tqdm(Pybel_object_dict[filename].items()):
        
        vapour_pressure_dict_BP[filename][species_name]['VP_Nannoolal']['BP_Nannoolal']=vapour_pressures.nannoolal(Pybel_object, temperature, boiling_point_dict[filename][species_name]['nannoolal'])
        vapour_pressure_dict_BP[filename][species_name]['VP_Nannoolal']['BP_Stein_Brown']=vapour_pressures.nannoolal(Pybel_object, temperature, boiling_point_dict[filename][species_name]['stein_and_brown'])
        vapour_pressure_dict_BP[filename][species_name]['VP_Nannoolal']['BP_Joback_Reid']=vapour_pressures.nannoolal(Pybel_object, temperature, boiling_point_dict[filename][species_name]['joback_and_reid'])
        vapour_pressure_dict_BP[filename][species_name]['VP_Myrdal_Yalkowsky']['BP_Nannoolal']=vapour_pressures.myrdal_and_yalkowsky(Pybel_object, temperature, boiling_point_dict[filename][species_name]['nannoolal'])
        vapour_pressure_dict_BP[filename][species_name]['VP_Myrdal_Yalkowsky']['BP_Stein_Brown']=vapour_pressures.myrdal_and_yalkowsky(Pybel_object, temperature, boiling_point_dict[filename][species_name]['stein_and_brown'])
        vapour_pressure_dict_BP[filename][species_name]['VP_Myrdal_Yalkowsky']['BP_Joback_Reid']=vapour_pressures.myrdal_and_yalkowsky(Pybel_object, temperature, boiling_point_dict[filename][species_name]['joback_and_reid'])
        vapour_pressure_dict[filename][species_name]['EVAP']=vapour_pressures.evaporation(Pybel_object, temperature)
        vapour_pressure_dict[filename][species_name]['EVAP2']=vapour_pressures.evaporation2(Pybel_object, temperature)
        vapour_pressure_dict[filename][species_name]['SIMPOL']=vapour_pressures.simpol(Pybel_object, temperature)

    
    # 2c) Subooled liquid density [(g/cc) at a specific temperature]
    # Some require a pre-calculated critical properties and thus boiling points which are fed into the function directly

    print('Liquid density, file:',filename)

    for species_name, Pybel_object in tqdm.tqdm(Pybel_object_dict[filename].items()):

        boiling_point=boiling_point_dict[filename][species_name]['nannoolal'] #Selecting this as an example

        density[filename][species_name]['Girolami']=np.real((liquid_densities.girolami(Pybel_object)))
    
        density_CP[filename][species_name]['Schroeder']['CP_Nannoolal']=np.real(liquid_densities.schroeder(Pybel_object, temperature, critical_properties.nannoolal(Pybel_object, boiling_point)))
        density_CP[filename][species_name]['Schroeder']['CP_Joback_Reid']=np.real(liquid_densities.schroeder(Pybel_object, temperature, critical_properties.joback_and_reid(Pybel_object, boiling_point)))
        density_CP[filename][species_name]['Le_Bas']['CP_Nannoolal']=np.real(liquid_densities.le_bas(Pybel_object, temperature, critical_properties.nannoolal(Pybel_object, boiling_point)))
        density_CP[filename][species_name]['Le_Bas']['CP_Joback_Reid']=np.real(liquid_densities.le_bas(Pybel_object, temperature, critical_properties.joback_and_reid(Pybel_object, boiling_point)))
        density_CP[filename][species_name]['Tyn_Calus']['CP_Nannoolal']=np.real(liquid_densities.tyn_and_calus(Pybel_object, temperature, critical_properties.nannoolal(Pybel_object, boiling_point)))
        density_CP[filename][species_name]['Tyn_Calus']['CP_Joback_Reid']=np.real(liquid_densities.tyn_and_calus(Pybel_object, temperature, critical_properties.joback_and_reid(Pybel_object, boiling_point)))

with open('density_test.json', 'w') as fp:
    json.dump(density, fp)

with open('density_cp_test.json', 'w') as fp:
    json.dump(density_CP, fp)

with open('bp_test.json', 'w') as fp:
    json.dump(boiling_point_dict, fp)

with open('vp_bp_test.json', 'w') as fp:
    json.dump(vapour_pressure_dict_BP, fp)

with open('vp_test.json', 'w') as fp:
    json.dump(vapour_pressure_dict, fp)