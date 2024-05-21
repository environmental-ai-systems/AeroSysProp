import pytest

##########################################################################################
# This file runs a series of tests based on stored outputs from existing methods.
#
#
# Current tests included:
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
import pdb

# Add the aerosysprop path
parent_directory = os.path.abspath('../..')
sys.path.append(parent_directory)

from aerosysprop.methods import boiling_points
from aerosysprop.methods import vapour_pressures
from aerosysprop.methods import critical_properties
from aerosysprop.methods import liquid_densities
