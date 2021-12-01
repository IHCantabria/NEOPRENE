'''
Library containing classes for computing statistics from time series.
The statistics are then used to calibrate the model parameters

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''


from NEOPRENE.NSRP.libs_NSRP import *
import numpy as np
import pandas as pd


def statistics_from_file(statistics_name, Seasonality_str, Seasonality, temporal_resolution, file):
    statististics_dataframe=pd.DataFrame(index=statistics_name,columns=Seasonality_str)
            
    statistics_external = pd.read_csv(file,index_col=0)
    statistics_external.columns = Seasonality
            
    for pri, prii in enumerate(Seasonality):
        statististics_values_real=statistics_external.iloc[:,pri].values
        statististics_dataframe.loc[:,str(prii)]=statististics_values_real
    return statististics_dataframe
    
def statistics_from_serie(statistics_name, Seasonality_str, Seasonality, temporal_resolution, time_series):
    statististics_dataframe=pd.DataFrame(index=statistics_name,columns=Seasonality_str)
    Datos_ = pd.DataFrame(index=time_series.index,columns = ['Rain'])
    Datos_.loc[:,'Rain'] = time_series.iloc[:,0].values
        
    for pri, prii in enumerate(Seasonality):
        Datos=Datos_.copy(); Datos[Datos['Rain']<0]=np.nan
             
        if len(Seasonality)==12:
            ## We select only the dates (months) for which I am going to calculate the statistics to be adjusted later.
            Datos['Estacionalidad']=Datos.loc[:,'Rain']*np.nan
            pos=np.where(Datos.index.month == prii); pos=pos[0]
            Datos['Estacionalidad'][pos]=Datos['Rain'][pos]
            Pluvio_GS = Datos['Estacionalidad'][Datos['Estacionalidad']>=0]
            Datos=Pluvio_GS
            
        else:
            ## We select only the dates (months) for which I am going to calculate the statistics to be adjusted later.
            Datos['Estacionalidad'] = Datos.loc[:,'Rain']*np.nan
            for i, ii in enumerate(prii):
                pos=np.where(Datos.index.month == ii); pos=pos[0]
                Datos['Estacionalidad'][pos]=Datos['Rain'][pos]
            Pluvio_GS = Datos['Estacionalidad'][Datos['Estacionalidad']>=0]
            Datos=Pluvio_GS
        
        
            ## We calculate the defined statistics to be adjusted and I enter them in a dataframe.
        statististics_values = calculate_statistics(Datos,statistics_name, temporal_resolution)
        
        statististics_dataframe.loc[:,str(prii)]=statististics_values
    return statististics_dataframe


class Statistics (object):
    """ This function allows the computation of the point model's spline statistics for the calibration from different data sources. 
    The input data can be a file with the statistics or a time series from which the statistics are calculated.
    
    Parameters
    ----------
    serie : Dataframe or None. 
        Precipitation time series with a given resolution. In the case where the series is not available and only the statistics are available, it is necessary to define this parameter as None
        
    file : string or None
        File containing the database containing the statistics previously defined in the input file containing the parameters for calibration.
        
    Outputs
    -------
    Dataframe with calculated statistics. 
    """
    def __init__(self, hiperparams, time_series=None, file=None):
    
        statistics_name     = hiperparams.statistics_name
        Seasonality_str     = hiperparams.Seasonality_str
        Seasonality         = hiperparams.Seasonality
        temporal_resolution = hiperparams.temporal_resolution
            
        if type(time_series)!=type(None) and file != None:
            raise Exception ('It is not possible to enter two types of data')
        if type(time_series)!=type(None):
            self.statististics_dataframe = statistics_from_serie(statistics_name,Seasonality_str,Seasonality,temporal_resolution, time_series)
        elif file != None:
            self.statististics_dataframe = statistics_from_file(statistics_name,Seasonality_str,Seasonality,temporal_resolution, file)
        else:
            raise Exception