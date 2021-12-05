
import sys
sys.path.insert(0, '/home/javi/source/doctorado/Point_Processes_Rain_models09102017/scikit-extremes')##importo la carpeta de extremos
#import skextremes2 as ske
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import datetime
from datetime import timedelta

def allmonths(dataframe):
    """Paso el dataframe a 12 columnas con 12 meses"""
    if len(dataframe.columns)==12: 
        dataframe_meses=dataframe

    else: 
        dataframe_meses=pd.DataFrame(index=dataframe.index)
        for i in dataframe.columns:
            if np.size(i)==1: dataframe_meses[i]=dataframe[i]
            else: 
                for ii in i: dataframe_meses[ii]=dataframe[i]

    dataframe=pd.DataFrame(index=dataframe_meses.index)
    for i in range(1, 13):
        dataframe[i]=dataframe_meses[i]
    return dataframe

def datetime2matlabdn(dt):
    mdn = dt + timedelta(days = 366)
    frac_seconds = (dt-datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds\
    /(24.0 * 60.0 * 60.0)
    frac_microseconds = dt.microsecond / (24.0 * 60.0 * 60.0 * 1000000.0)
    return mdn.toordinal() + frac_seconds + frac_microseconds



def datetime2matlabdnJavi(pdt, dt):
    '''
    Funcion que passa de un vectors de pandas a una lista de tiempo en formato ordinal
    pdt=pd.date_range
    dt = resolución days/hours/minutes/seconds
    '''
    
    if dt=='days':
        dt_ordinal=1
    elif dt=='hours':
        dt_ordinal=1/24
    elif dt=='minutes':
        dt_ordinal=1/(24*60)
    elif dt=='seconds':
        dt_ordinal=1/(24*60*60)
        
    d1 = datetime.datetime(pdt.year[0], pdt.month[0], pdt.day[0])# start date
    tt_ordinal_aux=list()
    
    dtt=0
    for i in range(len(pdt)):
        if i ==0:
            #tt_ordinal_origen=(d1).toordinal()
            tt_ordinal_origen=datetime2matlabdn(d1)
            tt_ordinal_aux.append(tt_ordinal_origen)
        else:
            tt_ordinal_aux.append(tt_ordinal_origen+dtt)
            
        dtt=dtt+dt_ordinal
        
    return tt_ordinal_aux

def RETURN_PERIODS(Dataframe, n_dias, tipo, ci_met, T):

    Dataframe_y = Dataframe.resample('A', how='max') #, how='max'
    x = Dataframe_y.values
    #########Selecciono aquellos datos donde hay mas de 292 dias de datos
    posi_mayores_80=list()
    for k, kk in enumerate(np.unique(Dataframe.index.year)):
        enc = Dataframe.index.year==kk; enc=np.sum(enc);
        if enc > n_dias:#292: #dias
            posi_mayores_80.append(k)
    x = x[posi_mayores_80]
    x = x[x>0]

    #Empiezo el ajuste	
    if tipo=='Gumbel':
        model = ske.models.classic.Gumbel(x, fit_method = 'mle', ci = 0.05, \
                    return_periods=T, ci_method = ci_met)
    elif tipo=='GEV':
        model = ske.models.classic.GEV(x, fit_method = 'mle', ci = 0.05, \
                    return_periods=T, ci_method = ci_met)
        
       

def func(x, a, b, c):
    """Funcion a la que voy a ajustar los valores de correlación espacial"""

    return a * np.exp(-b * x) + c

def distancia_f(x1, y1, x2, y2):
    dist=((x1-x2)**2 + (y1-y2)**2)**0.5
    return dist
#def haversine(lon1, lat1, lon2, lat2):
#    from math import radians, cos, sin, asin, sqrt
#    """
#    Calculate the great circle distance between two points 
#    on the earth (specified in decimal degrees)
#    """
#    # convert decimal degrees to radians 
#    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
#    # haversine formula 
#    dlon = lon2 - lon1 
#    dlat = lat2 - lat1 
#    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
#    c = 2 * asin(sqrt(a)) 
#    km = 6367 * c
#    return km

    model.params
    model.plot_summary()
    RETURN_PERIODS.return_periods=model.return_periods
    RETURN_PERIODS.return_values=model.return_values
    plt.show()

def IDW_f(xx, yy, zz, x, y, betha):
    ##xx, yy --> Coordenadas de las estaciones
    ##zz --> variable que quieres interpolar
    ##x, y --> punto donde quieres obtener el resultado
    ## betha--> exponente de peso. Cuanto mayor, más peso se le asigna a los puntos cercanos
    dist=list()
    for i in range(len(xx)):
        dist.append(haversine(xx[i], yy[i], x, y))
    dist=np.array(dist)
    weights=(dist**(-betha))/np.sum(dist**(-betha))
    
    return np.sum(zz*weights)
 

	

