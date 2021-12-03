'''
Library containing the main functions for the execution of the model.

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''

import pandas as pd

class inputs_simulation (object):
    def __init__(self,inputs_sim):
            
        if type(inputs_sim) != str and type(inputs_sim)!= list:
            self.params_cal = [inputs_sim.Fitted_parameters,inputs_sim.Dataframe_xi_months]
        elif type(inputs_sim) == str or type(inputs_sim) == list:
            self.params_cal = [pd.read_csv(inputs_sim[0],index_col=0), pd.read_csv(inputs_sim[1],index_col=0)]
        else:
            raise Exception