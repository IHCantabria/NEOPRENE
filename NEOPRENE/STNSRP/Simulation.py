'''
Library containing classes for simulating time series from the model.

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''



from NEOPRENE.NSRP.MathematicalPropertiesNSRP import *
from NEOPRENE.NSRP.utils import *
from NEOPRENE.NSRP.libs_NSRP import *
from NEOPRENE.NSRP.inputs_simulation import *
from NEOPRENE.NSRP.outputs_simulation import *
import time

class Simulation(object):
    def __init__(self,hiperparams):
        self.hiperparams = hiperparams
        
    def __call__(self, params_cal, Input_Series, Input_Attr):
   
        statististics_sim_df=pd.DataFrame(index=self.hiperparams.statistics_name,columns=self.hiperparams.Seasonality_str)
        
        Df_params = inputs_simulation(params_cal).params_cal
        
        Dataframe_xi_months = XI_MONTHS(Input_Series, CAL1, CAL)
        
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
        
        [Df_sim_join_hour_, Df_sim_join_day_]=\
            STNSRP_simulation(Df_params, XX, YY, self.hiperparams.year_ini, self.hiperparams.year_fin, self.hiperparams.temporal_resolution, self.hiperparams.process,
                             self.hiperparams.coordinates,self.hiperparams.storm_radius, self.hiperparams.Seasonality)
        
        ## Real, adjusted and simulated statistical calculations.
        

        if self.hiperparams.temporal_resolution == 'd':
            Data=Df_sim_join_day_.copy()
        elif self.hiperparams.temporal_resolution == 'h':
            Data=Df_sim_join_hour_.copy()
        Data[Data['Rain']<0]=np.nan
        Data[Data['Rain']<0]=np.nan
        
        
        statististics_sim_df = Statistics(self.hiperparams, time_series = Data, attributes = Input_Attr)

        results = outputs_simulation(Df_sim_join_day_,Df_sim_join_hour_,statististics_sim_df, self.hiperparams.temporal_resolution)

        return results
