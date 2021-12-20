'''
Library containing classes for keeping the simulation output.

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''

import pandas as pd

class outputs_simulation (object):
    def __init__(self,df_sim_join_daily,df_sim_join_hour,statististics_sim_df,crosscorr_sim_df,crosscorr_dataframe_dist, temporal_resolution, Seasonality):
        self.Daily_Simulation          = df_sim_join_daily
        self.Hourly_Simulation         = df_sim_join_hour
        self.statististics_Simulated   = statististics_sim_df
        self.crosscorr_Simulated       = crosscorr_sim_df
        self.crosscorr_Simulated_dist  = crosscorr_dataframe_dist
        self.temporal_resolution       = temporal_resolution
        self.Seasonality               = Seasonality

    def save_files(self,path_output_files):
        
        list_crosscorr = list()
        for j in self.crosscorr_Simulated.keys():
            if 'crosscorr_' in j:
                list_crosscorr.append(j)
        
        self.Daily_Simulation.to_csv(path_output_files+'Time_serie_daily_simulated.csv')
        self.Hourly_Simulation.to_csv(path_output_files+'Time_serie_hourly_simulated.csv')
        
        if self.temporal_resolution == 'd':
            for i in self.Seasonality:
                self.statististics_Simulated[i].to_csv(path_output_files+'statistic_daily_simulated_'+ str(i)+ '.csv')
                if len(list_crosscorr)!=0:
                    for c in list_crosscorr:
                        self.crosscorr_Simulated[c][i].to_csv(path_output_files+c+'_daily_simulated_'+ str(i)+ '.csv')
            
        elif self.temporal_resolution == 'h':
            for i in self.Seasonality:
                self.statististics_Simulated[i].to_csv(path_output_files+'statistic_hourly_simulated_'+ str(i)+ '.csv')
                if len(list_crosscorr)!=0:
                    for c in list_crosscorr:
                        self.crosscorr_Simulated[c][i].to_csv(path_output_files+c+'_hourly_simulated_'+ str(i)+ '.csv')
