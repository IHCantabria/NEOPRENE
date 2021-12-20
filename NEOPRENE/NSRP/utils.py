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
    """We move the dataframe to 12 columns with 12 months"""
    
    dataframe_meses=pd.DataFrame(index=dataframe.index, columns=np.arange(1,13))
    if len(dataframe.columns)==12:
        dataframe_meses.loc[:,:] = dataframe.values
    else:
        for i,ii in enumerate(dataframe.columns):
            for j in list(map(int, ii[1:-1].split(','))):
                dataframe_meses.loc[:,j] = dataframe.loc[:,ii]

   
    return dataframe_meses


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
def haversine(lon1, lat1, lon2, lat2):
    from math import radians, cos, sin, asin, sqrt
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    return km

    model.params
    model.plot_summary()
    RETURN_PERIODS.return_periods=model.return_periods
    RETURN_PERIODS.return_values=model.return_values
    plt.show()

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

def compare_statistics(CAL, SIM, frecuency):
    
    stats_obs=allmonths(CAL.statististics_Real)
    stats_fit=allmonths(CAL.statististics_Fit)
    stats_sim=allmonths(SIM.statististics_Simulated)
    
    N = len(stats_obs.index)
    RK = int(np.ceil(np.sqrt(len(stats_obs.index)))); 
    CK = int(math.ceil(N / RK))
        
    delete =  RK*CK-len(stats_obs.index)
           
    fig, axes = plt.subplots(RK, CK, figsize=(15, 15))
    
    axes = axes.ravel()
    for i, ii in enumerate(stats_obs.index):
        if 'mean' in ii:
            name_statistics = r'$\mu$'
        elif 'var' in ii:
            name_statistics = r'$\sigma_{'+ii.split('_')[-1]+frecuency+'}$'
        elif 'autocorr' in ii:
            name_statistics = r'$ACF-lag1_{'+ii.split('_')[-2]+frecuency+'}$'
        elif 'fih' in ii:
            name_statistics = r'$\phi_{'+ii.split('_')[-1]+frecuency+'}$'
        elif 'fiWW' in ii:
            name_statistics = r'$\phi^{WW}_{'+ii.split('_')[-1]+frecuency+'}$'
        elif 'fiDD' in ii:
            name_statistics = r'$\phi^{DD}_{'+ii.split('_')[-1]+frecuency+'}$'
        elif 'M3' in ii:
            name_statistics = r'$\overline{\mu}_{3_'+ii.split('_')[-1]+frecuency+'}$'
        
        Data_sta=pd.DataFrame(index=np.arange(1, 13))
        Obs=list(); Fit=list(); Sim=list();
        for m in np.arange(1, 13):
            Obs.append(stats_obs[m].loc[ii])
            Fit.append(stats_fit[m][ii])
            Sim.append(stats_sim[m].loc[ii])
            
        Data_sta.index = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

        Data_sta['Obs']=Obs
        Data_sta['Fit']=Fit
        Data_sta['Sim']=Sim
        #Ploteo
        ax=axes[i]
        ax.set_title(name_statistics)
        ax.grid(True)
        legnd=['Obs', 'Fit', 'Sim']
        pp=Data_sta['Obs'].plot(style='k--', lw=2,  ax=axes[i])
        Data_sta['Fit'].plot(style='bs', ax=axes[i])
        Data_sta['Sim'].plot(style='r^', ax=axes[i])
        #pp.legend(legnd, loc='best', )
        if ii[0:2]==('fi' or 'au'):
            ax.set_ylim(0, 1)
            
        ax.set_xticks(np.arange(0,12))
        ax.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],rotation = 45)
            
        ax.grid()
            
        del name_statistics     
    if delete!=0:
        for i in range(1,delete+1):
            axes[-i].remove()
        
    lines = []
    labels = []
    for i,ax in enumerate(fig.axes):
        axLine, axLabel = ax.get_legend_handles_labels()
        lines.extend(axLine)
        labels.extend(axLabel)
    fig.legend(lines[:3], labels[:3],           
               loc = 8,ncol=3,fontsize=15)
    fig.tight_layout(pad=3.7)
    
    return fig

def exceedence_probability(Serie_Observed, Serie_Simulated, temporal_resolution):
    # Observed
    aux = Serie_Observed.values
    aux = aux[aux>=0]
    s_obs = np.sort(aux, axis=None)[::-1]
    exc_obs = np.arange(1.,len(s_obs)+1) / len(s_obs)

    # Simulated
    if temporal_resolution == 'd':
        units = 'mm/day'
    elif temporal_resolution == 'h':
        units = 'mm/hour'

    aux = Serie_Simulated.values; aux=aux[aux>=0]
    s_sim = np.sort(aux, axis=None)[::-1]
    exc_sim = np.arange(1.,len(s_sim)+1) / len(s_sim)
    
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(exc_obs, s_obs, '.k', markersize=8,  label='Obs')
    ax.plot(exc_sim, s_sim, 'r',  label='Sim')
    ax.set_xlabel("Exceedence Probability", fontsize=13)
    ax.set_ylabel(units, fontsize=13)
    ax.set_xscale('log')
    ax.grid(True)
    ax.tick_params(labelsize=13)
    ax.legend(fontsize = 15)
    
    return fig

def disaggregate_rainfall(x_series, y_series):
    """
    Dissagregation function from: A spatial–temporal point process model of rainfall
    for the Thames catchment, UK (Cowpertwait 2005). Eq:15
    
    Inputs:
    x_series: Daily observed series (DataFrame)
    y_series: Hourly simulated series (DataFrame)
    
    Outputs
    results: x_series disaggregated from daily-to-hourly
    """

    #y_series_daily = y_series.resample('D').agg(pd.Series.sum, min_count=1)
    #results=x_series.resample('h').agg(pd.Series.sum, min_count=1)*np.nan
    
    y_series_daily = y_series.resample('D').agg(pd.Series.sum, min_count=1)
    dti = pd.date_range(start=x_series.index[0], end=x_series.index[-1] + timedelta(hours=23), freq="H")
    #results=pd.DataFrame(index = dti, columns = ['Rain'])
    results=pd.DataFrame(index = dti, columns = y_series.columns)

    for n, date in enumerate(x_series.index[1:]):
        x_j=x_series.loc[date].values       

        if np.nansum(x_j)==0:#day with no rain
            posi_d=np.where((results.index.year==date.year) & (results.index.month==date.month) & (results.index.day==date.day))[0]
            results.iloc[posi_d]=0
        else:
            x_j_1=x_series.loc[x_series.index[n]].values
            y_j_AUX=y_series_daily.iloc[1:]
            y_j_1_AUX=y_series_daily.iloc[0:-1]

            poss=np.sum(y_j_AUX[y_j_AUX.columns[np.where(x_j>0)[0]]]>0, axis=1)==len(np.where(x_j>0)[0])#days where observed and simulated rain is equal

            y_j=y_j_AUX.loc[poss].copy()
            y_j_1=y_j_1_AUX.iloc[np.where(poss==True)[0]].copy()        

            suma=np.nansum(((x_j-y_j.values)**2)+((x_j_1-y_j_1.values)**2), axis=1)

            posi=np.argmin(suma)
            posi_day=y_j.index[posi]
            rr_day=y_j.values[posi,:]

            posi_d=np.where((results.index.year==date.year) & (results.index.month==date.month) & (results.index.day==date.day))[0]
            posi_h_s=np.where((y_series.index.year==posi_day.year) & (y_series.index.month==posi_day.month)\
                              & (y_series.index.day==posi_day.day))[0]

            results.iloc[posi_d]=y_series.iloc[posi_h_s].values
            
    return results

def figure_disaggregation(hourly_disaggregation,daily_disaggregation,real_series,yearmin,yearmax):
    import matplotlib as mpl
    from matplotlib.lines import Line2D
    import matplotlib.patches as mpatches
    
    t1 = str(yearmin) + '-02-01'; t2 = str(yearmax) + '-02-28'
    
    hourly_disaggregation_c = hourly_disaggregation.loc[t1:t2].copy()
    daily_disaggregation_c  = daily_disaggregation.loc[t1:t2].copy()
    x_series_c              = real_series.loc[t1:t2].copy()
    
    hourly_disaggregation_tab = pd.DataFrame(index=hourly_disaggregation_c.index,columns=['mm','day','Hour','Par_Impar'])

    hourly_disaggregation_tab.iloc[:,0] = hourly_disaggregation_c.values
    hourly_disaggregation_tab.iloc[:,1] = hourly_disaggregation_c.index.day
    hourly_disaggregation_tab.iloc[:,2] = hourly_disaggregation_c.index.hour
    hourly_disaggregation_tab.iloc[:,3] = (hourly_disaggregation_tab.day.values % 2)
    
    fig, ax = plt.subplots(figsize=(20, 5))
    ax2 = ax.twinx()
   
    f    = x_series_c.iloc[:,0].plot(color = 'b', style= '--',  ax=ax,label = 'Obs. (daily)')
    h_im = (hourly_disaggregation_tab[(hourly_disaggregation_tab.Par_Impar==1)&((hourly_disaggregation_tab.mm>0.1))].mm.dropna()).plot(color = 'k', style = 's', markersize = 5, ax = ax2, label = 'Disagg. (hourly) day odd')
    h_p  = (hourly_disaggregation_tab[(hourly_disaggregation_tab.Par_Impar==0)&((hourly_disaggregation_tab.mm>0.1))].mm.dropna()).plot(color = 'k', style = 'x', markersize = 5, ax = ax2,label = 'Disagg. (hourly) day even')
    h_l  = ax2.plot(hourly_disaggregation_c.index,hourly_disaggregation_c.values, '--', c='grey', label = 'Disagg. (hourly)')

    l    = daily_disaggregation_c.iloc[:,0].plot(color = 'r', style = '^', markersize = 8, ax = ax, label = 'Disagg. (daily)')

    # ax.legend(['Obs. (daily)','Disagg. (daily)'], fontsize = 15)
    # ax2.legend(['Disagg. (hourly) day odd','Disagg. (hourly) day even','Disagg. (hourly)'], fontsize = 15, loc=2)
    ax.set_xlim(t1, t2)
    #ax2.set_xlim(t1, t2)
    ax.tick_params(axis = 'both', labelsize = 15)
    ax2.tick_params(axis = 'both', labelsize = 15)
    ax.set_yticks(np.arange(0,np.max(daily_disaggregation_c.values)+8,4))
    ax.set_ylim(0, np.max(daily_disaggregation_c.values)+8)
    ax2.set_ylim(0, (np.max(daily_disaggregation_c.values)+8)/2)
    ax.set_ylabel('mm/day', fontsize = 15);
    ax2.set_ylabel('mm/hour', fontsize = 15);
    ax.set_xticks(daily_disaggregation_c.index)
    ax.set_xticklabels(np.arange(1,len(daily_disaggregation_c)+1))
    ax.grid()
    ax2.grid()
    
    legend_elements = [Line2D([0], [0], label='Obs. (daily)', color='b',linestyle = '--'),
                   Line2D([], [], label='Disagg. (daily)', color='r',marker='^',linestyle='None'),
                   Line2D([0], [0], label='Disagg. (hourly) day odd', color='k',marker='s',linestyle = '--'),
                   Line2D([0], [0], label='Disagg. (hourly) day even', color='k',marker='x',linestyle = '--')]

    fig.legend(handles=legend_elements, loc = 8, ncol=5,fontsize=15)
    
    fig.tight_layout(pad=3.7)
    
    return fig
 

	

