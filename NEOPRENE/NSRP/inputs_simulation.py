'''
Library containing the main functions for the execution of the model.

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''

import pandas as pd

class inputs_simulation (object):
    def __init__(self,inputs_cal):
            
        if type(inputs_cal) != str:
            self.params_cal = inputs_cal.Fitted_parameters
        elif type(inputs_cal) == str:
            self.params_cal = pd.read_csv(inputs_cal,index_col=0)
        else:
            raise Exception