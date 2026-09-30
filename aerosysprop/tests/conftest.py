import os

# The suite runs at the default (full) guidance level whatever the shell has
# set. This must happen before the MCP server module is imported, since it
# fixes its level on import; tests that need bare set it themselves.
os.environ.pop('AEROSYSPROP_GUIDANCE', None)
