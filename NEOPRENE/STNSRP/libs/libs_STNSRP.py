'''
Librería que contine las funciones necesarias para simular el proceso espacio temporal de Neyman-Scott Rectangular Pulses (STNSRP)

	Autor: Javier Díez Sierra
	Fechas 06-02-2017
'''




'''Funciones del modelo NSRP'''

from MathematicalPropertiesSTNSRP import *
from libs_others import *
import PSOMdJ as pso
import numpy as np
import pandas as pd
import scipy as sp
import math as mt
import datetime
import scipy
from scipy.interpolate import interp1d
from datetime import date
from datetime import timedelta
from scipy import stats
from scipy.optimize import curve_fit

def calculate_statistics(Datos,statistics,temporal_resolution):
    """Calculo los estaísticos para una estacion determinada
    Los estadísticos que se pueden calcular son:
                statistics=[mean, var_h, autocorr_l_h, fih_h, fiWW_h, fiDD_h, M3_h]
                        (l=lag, h=aggregation levels)

    """
    
    statististics_values_real={}
    if temporal_resolution=='d':
        t='D'
    elif temporal_resolution=='h':
        t='h'

    #if np.sum(['mean' in i for i in statistics])>=1:
    #   statististics_values_real['mean']=np.nanmean(Datos)
    if np.sum(['var' in i for i in statistics])>=1:
        pos=np.where(['var' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",2)[1])
            aux=Datos.resample(str(h) + t , how='sum')
            statististics_values_real['var_'+str(h)]=sqrt(np.nanvar(aux))/np.nanmean(aux)
    if np.sum(['autocorr' in i for i in statistics])>=1:
        pos=np.where(['autocorr' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            l=int(statistics[ii].split("_",3)[1])
            h=int(statistics[ii].split("_",3)[2])
            aux=Datos.resample(str(h) + t , how='sum')
            Autocorrelation_aux=aux.autocorr(lag=l) 
            if np.size(Autocorrelation_aux)>1: Autocorrelation_aux=Autocorrelation_aux[0] 
            statististics_values_real['autocorr_'+str(l)+'_'+str(h)]=Autocorrelation_aux ##La correlación es igual a la covarianza entre la desviación tipica. 
            #yo calculo directamente la correlacion. (Aunque el el puntual lo estoy haciendo bien, en realidad debería de ajustar la covarianza)
    if np.sum(['fih' in i for i in statistics])>=1:
        pos=np.where(['fih' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",1)[1])
            statististics_values_real['fih_'+str(h)]=fi_h(Datos, h)
    if np.sum(['fiWW' in i for i in statistics])>=1:
        pos=np.where(['fiWW' in i for i in statistics]); pos=pos[0] 
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",1)[1])
            statististics_values_real['fiWW_'+str(h)]=fi_WW(Datos, h)
    if np.sum(['fiDD' in i for i in statistics])>=1:
        pos=np.where(['fiDD' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",1)[1])
            statististics_values_real['fiDD_'+str(h)]=fi_DD(Datos, h)	
    if np.sum(['M3' in i for i in statistics])>=1:
        pos=np.where(['M3' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",2)[1])
            aux=Datos.resample(str(h) + t , how='sum') 
            statististics_values_real['M3_'+str(h)]=(sp.stats.moment(aux, moment=3, nan_policy='omit'))/(np.nanvar(aux)**(3/2)) ##(aux o datos ¿?)

    return statististics_values_real

def statistics_per_station_and_stationality_f(Datos_, statistics, temporal_resolution, Seasonality):
    
    statistics_real_stationality_estations={}
    Datos=Datos_.copy()
    for i in Seasonality:
        statistics_stationality_estations_real=pd.DataFrame(index=statistics, columns=Datos.columns)
        for es in Datos.columns:
            ##Calculo los meses que no estan del 1 al 12. Pongo nans a los meses en los que no quiero calcular los
            ##estadísticos
            m_T=list()
            for n in np.arange(1, 13): 
                if not(np.sum(i==n)>=1): m_T.append(n)
            if len(m_T)==0: Datos_month=Datos[es].copy()
            elif len(m_T)==1: 
                Datos_month=Datos[es].copy()
                Datos_month[Datos_month.index.month==m_T]=np.nan
            else:
                Datos_month=Datos[es].copy(); 
                for m_t in m_T:
                    Datos_month[Datos_month.index.month==m_t]=np.nan

            Datos_month=Datos_month[Datos_month>=0]
            aux=calculate_statistics(Datos_month, statistics, temporal_resolution)

            #statistics_stationality_estations_real[es]=aux
            for ss in aux:
                statistics_stationality_estations_real[es][ss]=aux[ss]

        statistics_real_stationality_estations[i]=statistics_stationality_estations_real

    return statistics_real_stationality_estations



def cross_corr_stationality_f(DATOS, Seasonality, Estaciones, func, coordinates, cross_corr_h, temporal_resolution):
    
    cross_corr_stationality={}

    for i in Seasonality:
        Datos=DATOS.copy()
        cross_corr=pd.DataFrame()

        ##Calculo los meses que no estan del 1 al 12. Pongo nans a los meses en los que no quiero calcular los
        ##estadísticos
        m_T=list()
        for n in np.arange(1, 13): 
            if not(np.sum(i==n)>=1): m_T.append(n)
        if len(m_T)==0: Datos_month=Datos
        elif len(m_T)==1: 
            Datos_month=Datos
            Datos_month[Datos_month.index.month==m_T]=np.nan
        else:
            Datos_month=Datos; 
            for m_t in m_T:
                Datos_month[Datos_month.index.month==m_t]=np.nan
                
        h=int(cross_corr_h.split("_",3)[2])
        aux_month=Datos.resample(str(h) + temporal_resolution , how='sum')
        cross_corr_stationality_f.aa=aux_month

        cross_corr['dist'], cross_corr['cross_corr']=cross_correlation(Estaciones, aux_month, func, 11, coordinates)
        
        cross_corr_stationality[i]=cross_corr

    return cross_corr_stationality



def cross_correlation(Estaciones, Series, funcion, divisions, coordinates):
    Correlacion_espacial=pd.DataFrame(index=Series.columns, columns=Series.columns)
    Distancia=pd.DataFrame(index=Series.columns, columns=Series.columns)

    for i, ii in enumerate(Correlacion_espacial.index):
        x1=Estaciones[Estaciones['ID']==ii].X.values; x1=x1[0]
        y1=Estaciones[Estaciones['ID']==ii].Y.values; y1=y1[0]
        data1=Series[ii].values
        for j, jj, in enumerate(Correlacion_espacial.columns):
            x2=Estaciones[Estaciones['ID']==jj].X.values; x2=x2[0]
            y2=Estaciones[Estaciones['ID']==jj].Y.values; y2=y2[0]
            datos2=Series[jj].values
            if coordinates=='Geographical':
                Distancia[ii][jj]=haversine(x1, y1, x2, y2)
            elif coordinates=='UTM':
                Distancia[ii][jj]=distancia_f(x1, y1, x2, y2)
                Distancia=Distancia/1000#Lo paso a kilometros
                
            pos_no_nans=np.intersect1d(np.where(~np.isnan(Series[ii].values)), np.where(~np.isnan(Series[jj].values)))
            corr_s=scipy.stats.pearsonr(Series[ii].values[pos_no_nans], Series[jj].values[pos_no_nans])
            Correlacion_espacial[ii][jj]=corr_s[0]
    
    
     
    
    Correlacion_distancia=pd.DataFrame()
    Correlacion_distancia['Corr']=np.reshape(Correlacion_espacial.values, np.size(Correlacion_espacial.values), 1)
    Correlacion_distancia['dist']=np.reshape(Distancia.values, np.size(Distancia.values), 1)
    
    ##Quito los repetidos (solo me quedo con una parate de la diagonal)
    Correlacion_distancia=pd.DataFrame(np.vstack({tuple(row) for row in Correlacion_distancia.values}), columns=['Corr', 'dist'])
    ##Quito filas = nan
    Correlacion_distancia = Correlacion_distancia[np.isfinite(Correlacion_distancia['Corr'])]

    
    #fig= plt.figure(figsize=(20,8))
    #plt.plot(Correlacion_distancia['dist'], Correlacion_distancia['Corr'].values, '.r')
    #legnd = ['', 'Correlacion media']
    #plt.ylabel('correlation')
    #plt.xlabel('distance[km]')
    #plt.title('Spatial Correlation funtion')    
    
    ##Divido en 10 tramos y calculo la correlación espacial media en esa distancia
    tramos=np.linspace(0, Correlacion_distancia['dist'].max(), divisions)#divisions=11
    puntos_medios=(tramos[1:] + tramos[:-1]) / 2
    Correlacion_media=list()
    for i in range(len(puntos_medios)):
        pos1=np.where(Correlacion_distancia['dist']>=tramos[i]) 
        pos2=np.where(Correlacion_distancia['dist']<tramos[i+1])
        pos=np.intersect1d(pos1, pos2)
        Correlacion_media.append(np.mean(Correlacion_distancia['Corr'][pos]))    
    
    
    xdata = Correlacion_distancia['dist'].values
    ydata = Correlacion_distancia['Corr'].values
    
    cross_correlation.Correlacion_distancia=Correlacion_distancia
    
    popt, pcov = curve_fit(func, xdata, ydata)
    
    #plt.plot(xdata, ydata, 'b.')
    #plt.plot(xdata, func(xdata, popt[0], popt[1], popt[2]), 'r.')
    #plt.scatter(puntos_medios, func(puntos_medios, popt[0], popt[1], popt[2]), c='g', s=100)
    #sys.exit()
    
    puntos_medios
    correlacion_media=func(puntos_medios, popt[0], popt[1], popt[2])
    
    return puntos_medios, correlacion_media
