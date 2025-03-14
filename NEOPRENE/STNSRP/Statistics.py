'''
Library containing classes for computing statistics from time series.
The statistics are then used to calibrate the model parameters

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''

import sys
from NEOPRENE.STNSRP.libs_STNSRP import *
import numpy as np
import pandas as pd
import os


def statistics_from_file(statistics_name, Seasonality_str, Seasonality, temporal_resolution, files_folder):
    
    statistics_gauges={}
    
    crosscorr = list()
    for st in statistics_name:
        if 'crosscorr' in st:
            crosscorr.append(st)
            
    cross_corr_gauges={ccorr: {} for ccorr in crosscorr}
    
    
    files = os.listdir(files_folder)
    for sea in Seasonality:
        for file in files:
            if str(sea) in file and 'statististics' in file:
                stats = pd.read_csv(files_folder+file,index_col = 0)
                statistics_gauges[sea] = stats
                
    for cr in crosscorr:
        for sea in Seasonality:
            for file in files:
                if str(sea) in file and cr in file:
                    stats = pd.read_csv(files_folder+file,index_col = 0)
                    cross_corr_gauges[cr][sea] = stats
                    
    return statistics_gauges, cross_corr_gauges
    
def statistics_from_serie(statistics_name, Seasonality_str, Seasonality, temporal_resolution, time_series):
    statistics_name_no_cross = [sn for sn in statistics_name if 'crosscorr' not in sn]
    statistics_dataframe = pd.DataFrame(index=statistics_name_no_cross, columns=Seasonality_str)
    Datos_ = pd.DataFrame(index=time_series.index, columns=['Rain'])
    Datos_['Rain'] = pd.to_numeric(time_series.iloc[:, 0].values, errors='coerce')
    
    for prii in Seasonality:
        Datos = Datos_.copy()
        Datos[Datos['Rain'] < 0] = np.nan
             
        if len(Seasonality) == 12:
            pos_s = np.where(Datos_.index.month == prii)[0]
            pos_r = np.where(Datos_.values >= 0)[0]
            pos = np.intersect1d(pos_s, pos_r)
            Datos_time_period = pd.period_range('1800', periods=len(pos), freq=temporal_resolution)
            Datos = pd.DataFrame(Datos_['Rain'].iloc[pos].values, index=Datos_time_period)
        else:
            pos = []
            for ii in prii:
                pos_s = np.where(Datos_.index.month == ii)[0]
                pos_r = np.where(Datos_.values >= 0)[0]
                pos.append(np.intersect1d(pos_s, pos_r))
            pos = [item for sublist in pos for item in sublist]
            Datos_time_period = pd.period_range('1800', periods=len(pos), freq=temporal_resolution)
            Datos = pd.DataFrame(Datos_['Rain'].iloc[pos].values, index=Datos_time_period)
        
        statistics_values = calculate_statistics(Datos, statistics_name_no_cross, temporal_resolution)
        statistics_dataframe.loc[:, str(prii)] = statistics_values
    return statistics_dataframe

def statistics_from_several_series(statistics_name, Seasonality_str, Seasonality, temporal_resolution, time_series):
    statistics_gauges={}
    for ID in time_series.columns:
        time_serie = pd.DataFrame(index = time_series.index)
        time_serie['Rain'] = time_series[ID].values.astype(float)
        statistics_gauges[ID] = statistics_from_serie(statistics_name, Seasonality_str, Seasonality, temporal_resolution, time_serie)

    return statistics_gauges

def cross_corr(statistics_name, Seasonality, temporal_resolution, time_series, Attributes, func, coordinates):
    cross_corr_stationality = {}
    if any('cross' in i for i in statistics_name):
        pos = [i for i, val in enumerate(statistics_name) if 'cross' in val]
        for i in pos:
            cross_corr_s_aux, Distance_correlation = cross_corr_stationality_f(
                time_series, Seasonality, Attributes, func, coordinates, statistics_name[i], temporal_resolution
            )
            cross_corr_stationality[statistics_name[i]] = cross_corr_s_aux
    return cross_corr_stationality, Distance_correlation

def change_statistics_dic_order(statististics_dataframe, Seasonality):
    new_dic = {}
    keys = [k for k in statististics_dataframe.keys()]
    for sea in Seasonality:
        dtFr = pd.DataFrame(index = statististics_dataframe[keys[0]].index)
        for gauge in statististics_dataframe.keys():
            dtFr[gauge] = statististics_dataframe[gauge][str(sea)]#check
        new_dic[sea] = dtFr
    return new_dic


class Statistics (object):
    """ This function allows the computation of the point model's spline statistics and the spatal cross correlation for the calibration from different data sources. 
    The input data can be a file with the statistics or a time series from which the statistics are calculated.
    
    Parameters
    ----------
    time_series : Dataframe or None. 
        Precipitation time series with a given resolution. In the case where the series is not available and only the statistics are available, it is necessary to define this parameter as None

    attributes : Dataframe or None. 
        Attributes of the gauges (ID, Lon, Lat, elevation)
        
    file : string or None
        File containing the database containing the statistics previously defined in the input file containing the parameters for calibration.
        
    Outputs
    -------
    statististics_dataframe : Dataframe with calculated statistics. 
    crosscorr_dataframe : Dataframe with cross correlation.
    """
    def __init__(self, hiperparams, time_series=None, files_folder=None, attributes=None):
    
        statistics_name     = hiperparams.statistics_name
        Seasonality_str     = hiperparams.Seasonality_str
        Seasonality         = hiperparams.Seasonality
        temporal_resolution = hiperparams.temporal_resolution
        coordinates         = hiperparams.coordinates
            
        if type(time_series)!=type(None) and files_folder != None:
            raise Exception ('It is not possible to enter two types of data')
        if type(time_series)!=type(None) and type(attributes)!=type(None):
            if all(time_series.columns!=attributes.ID):
                raise Exception ('Lenght/name/order of time series and attributes does not match')
            else:
                #self.statististics_dataframe = statistics_from_serie(statistics_name,Seasonality_str,Seasonality,temporal_resolution, time_series)
                statististics_dataframe_dic = statistics_from_several_series(statistics_name,Seasonality_str,Seasonality,temporal_resolution, time_series)
                self.statistics_dataframe   = change_statistics_dic_order(statististics_dataframe_dic, Seasonality)
                self.crosscorr_dataframe, self.crosscorr_dist_dataframe     = cross_corr(statistics_name,Seasonality,temporal_resolution, time_series, attributes, func, coordinates)

        elif files_folder != None:
            [self.statististics_dataframe, self.crosscorr_dataframe] = statistics_from_file(statistics_name,Seasonality_str,Seasonality,temporal_resolution, files_folder)
        else:
            raise Exception