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
from openbabel import pybel
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
import pdb

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


print('Mapping species names to SMILES and Pybel objects to generate test data')

smiles_folder = '/SMILES'

SMILES_dict=collections.defaultdict(lambda: collections.defaultdict)
Pybel_object_dict=collections.defaultdict(lambda: collections.defaultdict)

#1) extract the data from each file and start generating a dictionary of Pybel objects
# The following is taken from PyBox
print('Mapping species names to SMILES and Pybel objects')

# In the following we try to account for the variable way in which SMILES strings are stored by application/project 
for filename in os.listdir(smiles_folder):

    if filename.endswith('.txt'):
        text=open(filename[:],'r')
        for line in text:
            input = line.split()
            try:
                # Keep a list of the information
                Compound_reference.append(input[0])
                smiles_array.append(input[1])
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
        dom = minidom.parse(SMILES_filename)
        for node in dom.getElementsByTagName('species'):
            species_name=node.attributes['species_name'].value
            if species_name not in ['OH','HO2','O3','hv']: 
                #Now check to see if the species is in your dictionary extracted from the mechanism
                if species_name in species_dict.values():
                    #pdb.set_trace()
                    #species_SMILES=node.attributes['smiles'].value
                    #Right, the following is an overkill way t extract the SMILES. If the XML was written correctly, and the SMILES were in quotation marks, the node. method would be all that is needed
                    try:
                        species_SMILES=re.search(r'<smiles>(.*?)</smiles>',node.getElementsByTagName('smiles').item(0).toxml()).group(1)
                        SMILES_dict[filename][species_name]=species_SMILES
                        Pybel_object=pybel.readstring('smi',species_SMILES)
                        Pybel_object_dict[filename][species_name]=Pybel_object
                    except:
                        print('No SMILES entry for species ',species_name)

# Save the SMILES dictionarites to a JSON file for later loading                        
with open('SMILES_master.json’, ‘w’) as fp:
          json.dump(SMILES_dict, fp)
# We need to save the Pybel objects in a different way
filehandler = open("Pybel_object_dict.obj","wb")
pickle.dump(Pybel_object_dict,filehandler)
filehandler.close()

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
    
