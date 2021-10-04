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
        
    def __call__(self, params_cal):
   
        statististics_sim_df=pd.DataFrame(index=self.hiperparams.statistics_name,columns=self.hiperparams.Seasonality_str)
        
        Df_params = inputs_simulation(params_cal).params_cal
        
        
        Df_params=allmonths(Df_params)
        print('')
        print('')
        print('#'*80)
        print('Synthetic simulation')
        print('')
        print('')
		## We start with synthetic simulation and statistical calculations.
        
        if self.hiperparams.process=='normal':
            [Df_sim_join_hour_, Df_sim_join_day_, Intensidad_cellss, Duracion_horas_cellss]=\
            NSRP_simulation(Df_params, self.hiperparams.year_ini, self.hiperparams.year_fin, self.hiperparams.temporal_resolution, self.hiperparams.process, self.hiperparams.Seasonality)
            print('Total cumulative rainfall - Analytical estimation = ' + "{:15.2f}".format(np.sum(Intensidad_cellss*Duracion_horas_cellss)))
            print('Total cumulative rainfall -             Simulated = ' + "{:15.2f}".format(np.sum(Df_sim_join_day_.values)))
            
        elif self.hiperparams.process=='storms':

            Df_params1=pd.DataFrame(index=['landa', 'ipsilon', 'eta', 'xi', 'betha'],columns=Df_params.columns)
            for i in Df_params.columns:
                Df_params1.loc[:,i]=Df_params.loc[:,i].values[0:5]

            Df_params2=pd.DataFrame(index=['landa', 'ipsilon', 'eta', 'xi', 'betha'],columns=Df_params.columns)
            for i in Df_params.columns:
                Df_params2.loc[:,i]=Df_params.loc[:,i].values[5:]

            [Df_sim_join_hour_1, Df_sim_join_day_1, Intensidad_cellss1, Duracion_horas_cellss1]=\
            NSRP_simulation(Df_params1, self.hiperparams.year_ini, self.hiperparams.year_fin, self.hiperparams.temporal_resolution, 'normal', self.hiperparams.Seasonality)

            [Df_sim_join_hour_2, Df_sim_join_day_2, Intensidad_cellss2, Duracion_horas_cellss2]=\
            NSRP_simulation(Df_params2, self.hiperparams.year_ini, self.hiperparams.year_fin, self.hiperparams.temporal_resolution, 'normal', self.hiperparams.Seasonality)

			##Combino los dos processos
            Df_sim_join_hour_=pd.DataFrame(index=Df_sim_join_hour_1.index)
            Df_sim_join_hour_['Rain']=Df_sim_join_hour_1.values+Df_sim_join_hour_2.values
            Df_sim_join_day_=pd.DataFrame(index=Df_sim_join_day_2.index)
            Df_sim_join_day_['Rain']=Df_sim_join_day_1.values+Df_sim_join_day_2.values

            print('Total cumulative rainfall - Analytical estimation - Storm 1 = ' + "{:15.2f}".format(np.sum(Intensidad_cellss1*Duracion_horas_cellss1))) 
            print('Total cumulative rainfall - Analytical estimation - Storm 2 = ' + "{:15.2f}".format(np.sum(Intensidad_cellss2*Duracion_horas_cellss2)))  

            print('Total cumulative rainfall - Analytical estimation = ' + \
            "{:15.2f}".format(np.sum(Intensidad_cellss1*Duracion_horas_cellss1)+np.sum(Intensidad_cellss2*Duracion_horas_cellss2)))
            print('Total cumulative rainfall -             Simulated = ' + "{:15.2f}".format(np.sum(Df_sim_join_day_.values)))    


        #statististics_dataframe=allmonths(statististics_dataframe)
        #print('')
        #print('')
        #print('#'*80)
        #print('Results Validation')
        #print('')
        #print('')
        
        ## Real, adjusted and simulated statistical calculations.

        for pri, prii in enumerate(self.hiperparams.Seasonality):
            if self.hiperparams.temporal_resolution == 'd':
                Data=Df_sim_join_day_.copy()
            elif self.hiperparams.temporal_resolution == 'h':
                Data=Df_sim_join_hour_.copy()
            Data[Data['Rain']<0]=np.nan
            Data[Data['Rain']<0]=np.nan

            if len(self.hiperparams.Seasonality)==12:
             
                Data['Seasonality']=Data['Rain']*np.nan
                pos=np.where(Data.index.month == prii); pos=pos[0]
                Data['Seasonality'].iloc[pos]=Data['Rain'].iloc[pos]
                Pluvio_GS = Data['Seasonality'][Data['Seasonality']>=0]
                Pluvio_GS[Pluvio_GS<0.001]=0
                Data=Pluvio_GS.astype(float)
                
            else:
  
                Data['Seasonality']=Data['Rain']*np.nan
                for i, ii in enumerate(prii):
                    pos=np.where(Data.index.month == ii); pos=pos[0]
                    Data['Seasonality'].iloc[pos]=Data['Rain'].iloc[pos]
                Pluvio_GS = Data['Seasonality'][Data['Seasonality']>=0]
                Pluvio_GS[Pluvio_GS<0.001]=0
                Data=Pluvio_GS.astype(float)



            statististics_values_synthetic = calculate_statistics(Data,self.hiperparams.statistics_name,self.hiperparams.temporal_resolution)
            statististics_sim_df.loc[:,str(prii)]=statististics_values_synthetic  

        #statististics_sim_df=allmonths(statististics_sim_df)
        
        results = outputs_simulation(Df_sim_join_day_,Df_sim_join_hour_,statististics_sim_df, self.hiperparams.temporal_resolution)

        return results
        
    # def save_files(self,results,path_output_files):
        
    #     results.Hourly_Simulation.to_csv(path_output_files+'Time_serie_hourly_simulated.csv')
    #     results.statististics_Simulated.to_csv(path_output_files+'statistic_hourly_simulated.csv')
    #     results.Daily_Simulation.to_csv(path_output_files+'Time_serie_day_simulated.csv')
