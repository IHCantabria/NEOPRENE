'''
Library containing classes for simulating time series from the model.
	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''



from NEOPRENE.STNSRP.MathematicalPropertiesSTNSRP import *
from NEOPRENE.STNSRP.utils import *
from NEOPRENE.STNSRP.libs_STNSRP import *
from NEOPRENE.STNSRP.inputs_simulation import *
from NEOPRENE.STNSRP.outputs_simulation import *
from NEOPRENE.STNSRP.Statistics import Statistics
import time

class Simulation(object):
    def __init__(self,hiperparams):
        self.hiperparams = hiperparams
        
    def __call__(self, params_cal, Input_Series, Input_Attr):
   
        statististics_sim_df=pd.DataFrame(index=self.hiperparams.statistics_name,columns=self.hiperparams.Seasonality_str)
    
        list_params = inputs_simulation(params_cal).params_cal
        
        Df_params = list_params[0]
        
        Df_params.columns = self.hiperparams.Seasonality
        
        Dataframe_xi_months = list_params[1]
        
        Df_params=allmonths(Df_params)
        
        XX = Input_Attr.X
        YY = Input_Attr.Y
        
        print('')
        print('')
        print('#'*80)
        print('Synthetic simulation')
        print('')
        print('')


		## We start with synthetic simulation and statistical calculations.

        if self.hiperparams.process=='normal':
   
            [Df_sim_join_day_,Df_sim_join_hour_]=\
                STNSRP_simulation(Df_params, Dataframe_xi_months, XX, YY, self.hiperparams.year_ini, self.hiperparams.year_fin, self.hiperparams.temporal_resolution, self.hiperparams.process,
                                 self.hiperparams.coordinates,self.hiperparams.storm_radius, self.hiperparams.Seasonality, Input_Attr.ID.values)

        elif self.hiperparams.process=='storms':

            if self.hiperparams.storm_radius == False:

                Df_params1=pd.DataFrame(index=['landa', 'ipsilon', 'eta', 'betha', 'fi_may'],columns=Df_params.columns)
                for i in Df_params.columns:
                    Df_params1.loc[:,i]=Df_params.loc[:,i].values[0:5]

                Df_params2=pd.DataFrame(index=['landa', 'ipsilon', 'eta', 'betha', 'fi_may'],columns=Df_params.columns)
                for i in Df_params.columns:
                    Df_params2.loc[:,i]=Df_params.loc[:,i].values[5:]

            elif self.hiperparams.storm_radius == True:

                Df_params1=pd.DataFrame(index=['landa', 'ipsilon', 'eta', 'betha', 'fi_may', 'fi_may_s'],columns=Df_params.columns)
                for i in Df_params.columns:
                    Df_params1.loc[:,i]=Df_params.loc[:,i].values[0:6]

                Df_params2=pd.DataFrame(index=['landa', 'ipsilon', 'eta', 'betha', 'fi_may', 'fi_may_s'],columns=Df_params.columns)
                for i in Df_params.columns:
                    Df_params2.loc[:,i]=Df_params.loc[:,i].values[6:]


            [Df_sim_join_day_1,Df_sim_join_hour_1]=\
            STNSRP_simulation(Df_params1, Dataframe_xi_months, XX, YY, self.hiperparams.year_ini, self.hiperparams.year_fin, self.hiperparams.temporal_resolution, 'normal',
                                 self.hiperparams.coordinates,self.hiperparams.storm_radius, self.hiperparams.Seasonality, Input_Attr.ID.values)

            [Df_sim_join_day_2,Df_sim_join_hour_2]=\
            STNSRP_simulation(Df_params2, Dataframe_xi_months, XX, YY, self.hiperparams.year_ini, self.hiperparams.year_fin, self.hiperparams.temporal_resolution, 'normal',
                                 self.hiperparams.coordinates,self.hiperparams.storm_radius, self.hiperparams.Seasonality, Input_Attr.ID.values)


            Df_sim_join_hour_=Df_sim_join_hour_1 + Df_sim_join_hour_2
            Df_sim_join_day_=Df_sim_join_day_1 + Df_sim_join_day_2 
        
        ## Real, adjusted and simulated statistical calculations.
        

        if self.hiperparams.temporal_resolution == 'd':
            Data=Df_sim_join_day_.copy()
        elif self.hiperparams.temporal_resolution == 'h':
            Data=Df_sim_join_hour_.copy()
            
        Data[Data<0.001]=0
    
        statististics_sim_df = Statistics(self.hiperparams, time_series = Data, attributes = Input_Attr)
        
        statistics_dataframe        = statististics_sim_df.statistics_dataframe.copy()
        crosscorr_dataframe         = statististics_sim_df.crosscorr_dataframe.copy()
        crosscorr_dataframe_dist    = statististics_sim_df.crosscorr_dist_dataframe.copy()

        results = outputs_simulation(Df_sim_join_day_,Df_sim_join_hour_,statistics_dataframe,crosscorr_dataframe,crosscorr_dataframe_dist, self.hiperparams.temporal_resolution,self.hiperparams.Seasonality)

        return results
