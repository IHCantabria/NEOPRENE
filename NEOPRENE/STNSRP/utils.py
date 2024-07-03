"""Library containing the complementary mathematical functions necessary for the operation of the code.

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
        
"""



import sys
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import datetime
from datetime import timedelta
import math

def allmonths(dataframe):
    """Transforma el dataframe a 12 filas con 12 meses."""
    expanded_rows = []

    # Expandir las tuplas en filas individuales
    for idx, row in dataframe.iterrows():
        # Convierte el índice de cadena a una tupla
        try:
            values = eval(idx)
        except:
            values = idx

        for val in values:
            # Agrega una nueva fila para cada valor en la tupla
            expanded_rows.append([val] + row.tolist())

    # Crear el nuevo DataFrame con las filas expandidas
    expanded_df = pd.DataFrame(expanded_rows, columns=["Index"] + list(dataframe.columns))

    # Establecer la columna 'Index' como índice
    expanded_df.set_index("Index", inplace=True)

    return expanded_df

def datetime2matlabdn(dt):
    mdn = dt + timedelta(days = 366)
    frac_seconds = (dt-datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds\
    /(24.0 * 60.0 * 60.0)
    frac_microseconds = dt.microsecond / (24.0 * 60.0 * 60.0 * 1000000.0)
    return mdn.toordinal() + frac_seconds + frac_microseconds



def datetime2matlabdnJavi(pdt, dt):
    '''
    Function that passes from a vector of pandas to a list of time in ordinal format.
    dt = resolutión days/hours/minutes/seconds
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

def func(x, a, b, c):
    """Function to which the spatial correlation values are fitted"""

    return a * np.exp(-b * x) + c

def distancia_f(x1, y1, x2, y2):
    dist=((x1-x2)**2 + (y1-y2)**2)**0.5
    return dist

def IDW_f(xx, yy, zz, x, y, betha):
    '''
    xx, yy --> Station coordinates
    zz --> variable to be interpolated
    x, y --> point where you want to get the result
    betha--> eweight exponent. The higher the exponent, the more weight is assigned to nearby points.
    '''
    dist=list()
    for i in range(len(xx)):
        dist.append(haversine(xx[i], yy[i], x, y))
    dist=np.array(dist)
    weights=(dist**(-betha))/np.sum(dist**(-betha))
    
    return np.sum(zz*weights)
 

	

