# vim: set et sw=4 sts=4 fileencoding=utf-8:
#
# Copyright (c) 2016 David Topping.
# All Rights Reserved.
# This file is part of AeroSysProp. It derives from UManSysProp
# (Topping et al., 2016, https://doi.org/10.5194/gmd-9-899-2016).
#
# AeroSysProp is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# AeroSysProp is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# AeroSysProp.  If not, see <http://www.gnu.org/licenses/>.


from . import data
from . import groups


def schroeder(compound):
    return groups.aggregate_matches(
        groups.schroeder(compound),
        data.SCHROEDER_DENSITY
        )


def le_bas(compound):
    return groups.aggregate_matches(
        groups.le_bas(compound),
        data.LE_BAS_DENSITY
        )


def tyn_and_calus(volume):
    return 0.285 * (volume ** 1.048)

