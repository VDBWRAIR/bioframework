__projectname__ = 'bioframework'
__release__ = '0.0.1'
__authors__ = "Tyghe Vallard, Michael Panciera"
__authoremails__ = 'vallardt@gmail.com, michael.panciera.work@gmail.com'
__description__ = "Basic building blocks for bio-pipelines"
__keywords__ = "pyjip, pipeline, bioinformatics"

from os.path import dirname, join, abspath
import os
JIP_PATH=join(dirname(dirname(abspath(__file__))), 'jip_modules')
