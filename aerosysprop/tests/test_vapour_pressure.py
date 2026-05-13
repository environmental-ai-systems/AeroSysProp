import pytest

##########################################################################################
# pytest will run all files of the form test_*.py or *_test.py in the current directory and its subdirectories.
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

step=0
filenames=[]
Compound_reference=[]
smiles_array=[]
property_array=[]

print('Testing properties against test [predicted] data')

# Load the relevant test data into dictionaries.
bp_test_f = open('bp_test.json') 
data_bp_test = json.load(bp_test_f) 

density_test_f = open('density_test.json') 
density_test = json.load(density_test_f) 

density_cp_test_f = open('density_cp_test.json') 
density_cp_test = json.load(density_cp_test_f) 

vp_bp_test_f = open('vp_bp_test.json') 
vp_bp_test = json.load(vp_bp_test_f) 

vp_test_f = open('vp_test.json') 
vp_test = json.load(vp_test_f) 


for filename in os.listdir(smiles_folder):
    
    # Now leed to loop through 'Pybel_object_dict[filename][species_name]=Pybel_object'
    print('Boiling points, file:',filename)

    for species_name, Pybel_object in tqdm.tqdm(Pybel_object_dict[filename].items()):

        assert(data_bp_test[filename][species_name]['joback_and_reid']==boiling_points.joback_and_reid(Pybel_object))
        assert(data_bp_test[filename][species_name]['stein_and_brown']==boiling_points.stein_and_brown(Pybel_object))
        assert(data_bp_test[filename][species_name]['nannoolal']==boiling_points.nannoolal(Pybel_object))
   
    # Now leed to loop through 'Pybel_object_dict[filename][species_name]=Pybel_object'
    print('Vapour pressure (psat), file:',filename)

    temperature=298.15

    for species_name, Pybel_object in tqdm.tqdm(Pybel_object_dict[filename].items()):
        
        assert(vp_bp_test[filename][species_name]['VP_Nannoolal']['BP_Nannoolal']==vapour_pressures.nannoolal(Pybel_object, temperature, data_bp_test[filename][species_name]['nannoolal']))
        assert(vp_bp_test[filename][species_name]['VP_Nannoolal']['BP_Stein_Brown']==vapour_pressures.nannoolal(Pybel_object, temperature, data_bp_test[filename][species_name]['stein_and_brown']))
        assert(vp_bp_test[filename][species_name]['VP_Nannoolal']['BP_Joback_Reid']==vapour_pressures.nannoolal(Pybel_object, temperature, data_bp_test[filename][species_name]['joback_and_reid']))
        assert(vp_bp_test[filename][species_name]['VP_Myrdal_Yalkowsky']['BP_Nannoolal']==vapour_pressures.myrdal_and_yalkowsky(Pybel_object, temperature, data_bp_test[filename][species_name]['nannoolal']))
        assert(vp_bp_test[filename][species_name]['VP_Myrdal_Yalkowsky']['BP_Stein_Brown']==vapour_pressures.myrdal_and_yalkowsky(Pybel_object, temperature, data_bp_test[filename][species_name]['stein_and_brown']))
        assert(vp_bp_test[filename][species_name]['VP_Myrdal_Yalkowsky']['BP_Joback_Reid']==vapour_pressures.myrdal_and_yalkowsky(Pybel_object, temperature, data_bp_test[filename][species_name]['joback_and_reid']))
        #pdb.set_trace()
        assert(vp_test[filename][species_name]['EVAP']==vapour_pressures.evaporation(Pybel_object, temperature))
        assert(vp_test[filename][species_name]['EVAP2']==vapour_pressures.evaporation2(Pybel_object, temperature))
        assert(vp_test[filename][species_name]['SIMPOL']==vapour_pressures.simpol(Pybel_object, temperature))

    
    # 2c) Subooled liquid density [(g/cc) at a specific temperature]
    # Some require a pre-calculated critical properties and thus boiling points which are fed into the function directly

    print('Liquid density, file:',filename)

    for species_name, Pybel_object in tqdm.tqdm(Pybel_object_dict[filename].items()):

        boiling_point=data_bp_test[filename][species_name]['nannoolal'] #Selecting this as an example

        assert(density_test[filename][species_name]['Girolami']==np.real((liquid_densities.girolami(Pybel_object))))
    
        assert(density_cp_test[filename][species_name]['Schroeder']['CP_Nannoolal']==np.real(liquid_densities.schroeder(Pybel_object, temperature, critical_properties.nannoolal(Pybel_object, boiling_point))))
        assert(density_cp_test[filename][species_name]['Schroeder']['CP_Joback_Reid']==np.real(liquid_densities.schroeder(Pybel_object, temperature, critical_properties.joback_and_reid(Pybel_object, boiling_point))))
        assert(density_cp_test[filename][species_name]['Le_Bas']['CP_Nannoolal']==np.real(liquid_densities.le_bas(Pybel_object, temperature, critical_properties.nannoolal(Pybel_object, boiling_point))))
        assert(density_cp_test[filename][species_name]['Le_Bas']['CP_Joback_Reid']==np.real(liquid_densities.le_bas(Pybel_object, temperature, critical_properties.joback_and_reid(Pybel_object, boiling_point))))
        assert(density_cp_test[filename][species_name]['Tyn_Calus']['CP_Nannoolal']==np.real(liquid_densities.tyn_and_calus(Pybel_object, temperature, critical_properties.nannoolal(Pybel_object, boiling_point))))
        assert(density_cp_test[filename][species_name]['Tyn_Calus']['CP_Joback_Reid']==np.real(liquid_densities.tyn_and_calus(Pybel_object, temperature, critical_properties.joback_and_reid(Pybel_object, boiling_point))))

