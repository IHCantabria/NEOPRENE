'''
Library containing classes for keeping the simulated model parameters.

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''


import pandas as pd

class parameters_calibration(object):
    def __init__(self, Dataframe_parametros, statististics_fit_dataframe, statististics_Real, Error_dataframe):
        self.Fitted_parameters  = Dataframe_parametros
        self.statististics_Fit  = statististics_fit_dataframe
        self.statististics_Real = statististics_Real
        self.Error_dataframe    = Error_dataframe

    def save_files(self, path_output_files): 
        self.Fitted_parameters.to_csv(path_output_files+'./Calibrated_parameters.csv')
        self.statististics_Real.to_csv(path_output_files+'./statististics_real.csv')
        self.statististics_Fit.to_csv(path_output_files+'./statististics_fit.csv')
            
        