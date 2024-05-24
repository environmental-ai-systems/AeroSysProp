import pytest

##########################################################################################
# pytest will run all files of the form test_*.py or *_test.py in the current directory and its subdirectories.
#
# Last modification 17/05/24

import glob
from openbabel import pybel
import collections
import sys
import os
import json
import pickle 
from xml.dom import minidom
import re
import tqdm
import json
import numpy as np

parent_directory = os.path.abspath('../..')
sys.path.append(parent_directory)

from aerosysprop.methods import boiling_points
from aerosysprop.methods import vapour_pressures
from aerosysprop.methods import critical_properties
from aerosysprop.methods import liquid_densities
import pdb

# 1) Generate dictionary of SMILES objects to pass to functions to generate results

##########################################################################################
# A) -  Individual property prediction
# We can select a number of SMILES string files. Some are sstored as XML files, others plain text
# The following loop works through the files stored in the SMILES folder 


print('Mapping species names to SMILES and Pybel objects to generate test data')

smiles_folder = os.getcwd()+'/SMILES'

SMILES_dict=collections.defaultdict(lambda: collections.defaultdict())
Pybel_object_dict=collections.defaultdict(lambda: collections.defaultdict())

#1) extract the data from each file and start generating a dictionary of Pybel objects
# The following is taken from PyBox
print('Mapping species names to SMILES and Pybel objects')

# In the following we try to account for the variable way in which SMILES strings are stored by application/project 
for filename in os.listdir(smiles_folder):

    if filename.endswith('.txt'):
        text=open(smiles_folder+'/'+filename,'r')
        for line in text:
            input = line.split()
            try:
                # Keep a list of the information
                #pdb.set_trace()
                #Compound_reference.append(input[0])
                #smiles_array.append(input[1])
                species_SMILES=input[1]
                # Now create Pybel objects which are used in all property predictive techniques
                # Note for text files, unlike XML files, its unlikely that we will also have a species name so using 'species_SMILES'
                SMILES_dict[filename][species_SMILES]=species_SMILES
                Pybel_object=pybel.readstring('smi',input[1])
                Pybel_object_dict[filename][species_SMILES]=Pybel_object
            except:
                print('Problem generating Pybel object from SMILES string')
    elif filename.endswith('.xml'):
        # parse an xml file by name - use this to store species name versus SMILES
        dom = minidom.parse(smiles_folder+'/'+filename)
        for node in dom.getElementsByTagName('species'):
            species_name=node.attributes['species_name'].value
            if species_name not in ['OH','HO2','O3','hv']: 
                #Now check to see if the species is in your dictionary extracted from the mechanism
                #Right, the following is an overkill way t extract the SMILES. If the XML was written correctly, and the SMILES were in quotation marks, the node. method would be all that is needed
                try:
                    #pdb.set_trace()
                    species_SMILES=re.search(r'<smiles>(.*?)</smiles>',node.getElementsByTagName('smiles').item(0).toxml()).group(1)
                    SMILES_dict[filename][species_name]=species_SMILES
                    Pybel_object=pybel.readstring('smi',species_SMILES)
                    Pybel_object_dict[filename][species_name]=Pybel_object
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

def test_bp():
    for filename in os.listdir(smiles_folder):
        # Now leed to loop through 'Pybel_object_dict[filename][species_name]=Pybel_object'
        print('Boiling points, file:',filename)
        for species_name, Pybel_object in tqdm.tqdm(Pybel_object_dict[filename].items()):
            assert(data_bp_test[filename][species_name]['joback_and_reid']==boiling_points.joback_and_reid(Pybel_object))
            assert(data_bp_test[filename][species_name]['stein_and_brown']==boiling_points.stein_and_brown(Pybel_object))
            assert(data_bp_test[filename][species_name]['nannoolal']==boiling_points.nannoolal(Pybel_object))
    return 

def test_vp():
    temperature=298.15
    for filename in os.listdir(smiles_folder):
        print('Psat, file:',filename)
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
    return

def test_dens():    
    # 2c) Subooled liquid density [(g/cc) at a specific temperature]
    # Some require a pre-calculated critical properties and thus boiling points which are fed into the function directly

    temperature=298.15
    for filename in os.listdir(smiles_folder):
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
    return 
