"""
1e. Using logging
===============================
LoopStructural has a number of levels of logging incorporated in the code to allow 
for recording and debugging the models.
The python logging module allows for 5 levels of messages to be returned to the user:
1. Debug messages
2. Info messages
3. Warning messages
4. Error messages
5. Critical messages

LoopStructural uses all of these logging levels to report the various aspects of the model
building process. 
Generally, the user only needs to be aware of the warning and error messages. 

By default the warning, error and critical messages are returned to the console and will appear to
the user. 
All messages except for debug are recorded to a file :code:`default-loop-structural-logfile.log`.

Lets have a look at the logging from the Claudius model. 
"""

from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer 
from LoopStructural.datasets import load_claudius #demo data 
from LoopStructural import log_to_file
import pandas as pd
import numpy as np

##################################################################################################
# Specify a log file
# ~~~~~~~~~~~~~~~~~~~~

log_to_file('logging_demo_log.log')

##################################################################################################
# Create model
# ~~~~~~~~~~~~~~~~~~~~
data, bb = load_claudius()
model = GeologicalModel(bb[0,:],bb[1,:])
model.set_model_data(data)

vals = [0,60,250,330,600]
strat_column = {'strati':{}}
for i in range(len(vals)-1):
    strat_column['strati']['unit_{}'.format(i)] = {'min':vals[i],'max':vals[i+1],'id':i}
model.set_stratigraphic_column(strat_column)
strati = model.create_and_add_foliation("strati",
                                           interpolatortype="FDI", # try changing this to 'PLI'
                                           nelements=1e4, # try changing between 1e3 and 5e4
                                           buffer=0.3,
                                           solver='pyamg',
                                           damp=True
                                          )
viewer = LavaVuModelViewer(model,background="white")
viewer.add_model_surfaces()
viewer.rotate([-85.18760681152344, 42.93233871459961, 0.8641873002052307])
viewer.display()
#################################################################################################
# Looking at the log file
# ~~~~~~~~~~~~~~~~~~~~~~~
# Here are the first 10 lines of the log file. 
# Most operations in loopstructural are recorded and this will allow you to identify whether 
# an operation is not occuring as you would expect.

from itertools import islice
# with open('logging_demo_log.log') as inf:
#     for line in islice(inf, 0, 11):
#         print(line)


#################################################################################################
# Logging to console
# ~~~~~~~~~~~~~~~~~~
# It is also possible to change the logging level for the console log.

from LoopStructural import log_to_console

log_to_console('info') 


from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer 
from LoopStructural.datasets import load_claudius #demo data 

import pandas as pd
import numpy as np


data, bb = load_claudius()
model = GeologicalModel(bb[0,:],bb[1,:])
model.set_model_data(data)

vals = [0,60,250,330,600]
strat_column = {'strati':{}}
for i in range(len(vals)-1):
    strat_column['strati']['unit_{}'.format(i)] = {'min':vals[i],'max':vals[i+1],'id':i}
model.set_stratigraphic_column(strat_column)
strati = model.create_and_add_foliation("strati",
                                           interpolatortype="FDI", # try changing this to 'PLI'
                                           nelements=1e4, # try changing between 1e3 and 5e4
                                           buffer=0.3,
                                           solver='pyamg',
                                           damp=True
                                          )
viewer = LavaVuModelViewer(model,background="white")
viewer.add_model_surfaces()
viewer.rotate([-85.18760681152344, 42.93233871459961, 0.8641873002052307])
viewer.display()