'''
Library containing classes for keeping the simulated model parameters.

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''


import pandas as pd

class parameters_calibration(object):
    def __init__(self, Dataframe_parametros, statististics_fit_dataframe, crosscorr_dataframe_fit, statististics_Real, crosscorr_Real, Error_dataframe, Dataframe_xi_months):
        self.Fitted_parameters  = Dataframe_parametros
        self.statististics_Fit  = statististics_fit_dataframe
        self.statististics_Real = statististics_Real
        self.crosscorr_Fit  = crosscorr_dataframe_fit
        self.crosscorr_Real = crosscorr_Real
        self.Error_dataframe    = Error_dataframe
        self.Dataframe_xi_months    = Dataframe_xi_months

    def save_files(self, path_output_files): 
        self.Fitted_parameters.to_csv(path_output_files   + 'Calibrated_parameters.csv')
        self.statististics_Fit.to_csv(path_output_files   + 'statististics_fit.csv')
        self.Dataframe_xi_months.to_csv(path_output_files + 'xi_months.csv')
        
        for i in self.Fitted_parameters.columns:
            self.statististics_Real[i].to_csv(path_output_files+'statististics_real_'+str(i)+'.csv')
            
            
        