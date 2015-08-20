#!/bin/bash -e

# F0401: Unable to import 'cosmos'
# W0612: unused variable in cmd string inside legos. They are used but pylint can't catch them.
# W0221: Arguments number differs from overridden 'cmd' method
# W0141: Used builtin function 'filter'  # This is a legit use of filter. Not sure why pylint complains about it
# W0110: map/filter on lambda could be replaced by comprehension.
# E1002: (super-on-old-class), BwaMem.__init__] Use of super on an old style class
# R0801: Similar lines
#ignored-modules=cosmos since pylint keep compalining about not being able to import cosmos
pylint --disable=F0401,E1120,W0612,W0221,W0141,E1002,W0110,R0801,W0613 ampsim/

# Disable F401: unused imports. Only really applies to __init__.py files, since all other unused imports
# Disable W391:  blank line at end of file but pylint requires a blank line at EOF!

flake8 --max-line-length=120 --ignore=F401 --ignore=W391 ampsim/