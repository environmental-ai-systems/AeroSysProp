# This file is part of AeroSysProp.
#
# AeroSysProp is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
"""
aerosysprop.guidance
====================
How much explanatory text the agent-facing layer serves.

``AEROSYSPROP_GUIDANCE`` selects the level:

* ``full`` (the default): tool descriptions, the agent guide resources, method
  card content in ``list_methods``, and the coverage statement, caveats and
  reporting requirements in batch summaries.
* ``bare``: tool names and arguments, method IDs, output columns and units.

The level changes only descriptive text. Values, diagnostics, provenance
labels, error messages and run folders are the same at every level, which is
what lets the agent evaluation attribute a difference in answers to the
documentation alone.
"""
from __future__ import annotations

import os

GUIDANCE_ENV = 'AEROSYSPROP_GUIDANCE'
LEVELS = ('full', 'bare')


def level() -> str:
    """The guidance level from the environment; unset or empty means full.

    Raises ValueError for an unknown level, so a typo cannot silently serve
    the wrong condition.
    """
    value = os.environ.get(GUIDANCE_ENV, '').strip().lower() or 'full'
    if value not in LEVELS:
        raise ValueError(f'{GUIDANCE_ENV}={value!r} is not a guidance level; '
                         f'use one of {", ".join(LEVELS)}.')
    return value


def is_bare() -> bool:
    return level() == 'bare'
