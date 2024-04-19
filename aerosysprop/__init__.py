# This file is part of AeroSysProp.
#
# AeroSysProp is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# AeroSysProp is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# umansysprop. If not, see <http://www.gnu.org/licenses/>.
"""An system for calculating the properties of individual organic
molecules and ensemble mixtures relevant to aerosol science"""
from __future__ import (
    unicode_literals,
    absolute_import,
    print_function,
    division,
    )
str = type('')
import sys

__project__ = 'aerosysprop'
__version__ = '0.1'
__keywords__ = ['aerosol', 'organic', 'chemistry','droplets','atmospheric']
__authors__ = 'David Topping and project team'
__lead_author_email__ = '[David Topping] david.topping@manchester.ac.uk'
__url__ = 'TBC'
__platforms__ = 'ALL'
__requires__ = ['openbabel',]
__extra_requires__ = {'server': ['openbabel'],
    'client': [],
    'doc': ['sphinx'],
    }
