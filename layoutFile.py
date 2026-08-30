"""Legacy layout compatibility constants.

The application UI moved to :mod:`ddmsoft_gui`. These names remain harmless
for old scripts that imported the layout module, but no GUI toolkit is loaded
from here.
"""

DEFAULTDATADIRECTORY = "~"
DEFAULTFONT = ("Sans Serif", 10)
TITLEFONT = ("Sans Serif", 12)
layout = None
