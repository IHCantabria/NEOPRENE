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

def compare_statistics(CAL, SIM, frecuency):

    dataframe_statististics_obs = pd.DataFrame(index = CAL.statististics_Fit.index)
    dataframe_statististics_fit = pd.DataFrame(index = CAL.statististics_Fit.index)
    dataframe_statististics_sim = pd.DataFrame(index = CAL.statististics_Fit.index)

    for season in CAL.statististics_Real:
        aux_season_obs = CAL.statististics_Real[season]
        aux_season_fit = CAL.statististics_Fit[str(season)]
        aux_season_sim = SIM.statististics_Simulated[season]
        dataframe_statististics_obs[season] = aux_season_obs.mean(axis=1)
        dataframe_statististics_fit[season] = aux_season_fit
        dataframe_statististics_sim[season] = aux_season_sim.mean(axis=1)
        
    stats_obs=allmonths(dataframe_statististics_obs)
    stats_fit=allmonths(dataframe_statististics_fit)
    stats_sim=allmonths(dataframe_statististics_sim)
    
     
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
        Data_sta.index = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
        Obs=list(); Fit=list(); Sim=list();
        for m in np.arange(1, 13):
            Obs.append(stats_obs[m].loc[ii])
            Fit.append(stats_fit[m][ii])
            Sim.append(stats_sim[m].loc[ii])

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
        #pp.legend(legnd, loc='best', fontsize=15)
        if ii[0:2]==('fi' or 'au'):
            ax.set_ylim(0, 1)
        del name_statistics
        
        ax.grid()
        ax.set_xticks(np.arange(0,12))
        ax.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],rotation = 45)
        
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
    
    ax.set_visible(False)
    
    return fig

def exceedence_probability(Serie_Observed, Serie_Simulated, temporal_resolution):
    # Observed
    aux=Serie_Observed.values; aux=aux[aux>=0]
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

    ax.plot(exc_obs*100, s_obs, '.k', markersize=8,  label='Obs')
    ax.plot(exc_sim*100, s_sim, 'r',  label='Sim')
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


def figure_correlation(CAL,SIM):

    import matplotlib as mpl
    from matplotlib.lines import Line2D
    import matplotlib.patches as mpatches

    RK=int(np.ceil(np.sqrt(len(SIM.Seasonality)))); 
    CK=int(np.ceil(np.sqrt(len(SIM.Seasonality)))); 
    
    figures = []
    names    = []

    for cross_sim in SIM.crosscorr_Simulated.keys():
        fig, axes = plt.subplots(RK, CK, figsize=(10, 10))
        for i, season in enumerate(SIM.Seasonality):
            obs_cross = CAL.crosscorr_Real[cross_sim][season]
            cal_cross = CAL.crosscorr_Fit[cross_sim][season]
            sim_cross = SIM.crosscorr_Simulated[cross_sim][season]

            tim = cross_sim.split('_')[-1]
            frq = SIM.temporal_resolution

            row = i // CK
            col = i % CK

            ax=axes[row,col]
            
            ax.plot(obs_cross['dist'].values, obs_cross['cross_corr'].values, 'b^-')
            ax.plot(cal_cross['dist'].values, cal_cross['cross_corr'].values, 'ks-')
            ax.plot(CAL.crosscorr_Real_Dist[season].dist, CAL.crosscorr_Real_Dist[season].Corr,'.b')
            ax.plot(sim_cross['dist'].values, sim_cross['cross_corr'].values, '-r')
            ax.plot(SIM.crosscorr_Simulated_dist[season].dist, SIM.crosscorr_Simulated_dist[season].Corr,'xr')
            ax.set_ylim(-0.1,1.1)
            ax.set_yticks(np.arange(0,1.1,0.2))
            #ax.set_yticklabels(np.arange(0,1.1,0.1))
            ax.set_title('Months '+str(season),fontsize = 15)
            ax.set_ylabel(tim+frq+' cross-correlation',fontsize = 15)
            ax.set_xlabel('distance (km)',fontsize = 15)
            ax.tick_params(axis = 'both', labelsize = 15)
            ax.grid()


        legend_elements = [Line2D([0], [0], label='Observed', color='b',linestyle = '-',marker='^'),
                           Line2D([0], [0], label='Calibrated', color='k',linestyle = '-',marker='s'),
                           Line2D([0], [0], label='Simulated', color='r',linestyle = '-'),
                           Line2D([], [], label='Simulated', color='r',marker='x',linestyle='None'),
                           Line2D([], [], label='Observed', color='b',marker='.',linestyle='None'),
                          ]

        fig.legend(handles=legend_elements, loc = 8, ncol=5,fontsize=15)

        fig.tight_layout(pad=3.7)
        
        figures.append(fig)
        names.append(cross_sim)
        
    return figures, names
 

	

