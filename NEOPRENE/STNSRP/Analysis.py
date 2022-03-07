'''
Library containing classes for calibrating model parameters.

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''

from NEOPRENE.STNSRP.utils import *


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
            #name_statistics = r'$\overline{\mu}_{3_'+ii.split('_')[-1]+frecuency+'}$'
            name_statistics = r'$\bar{\mu}_{3_{'+ii.split('_')[-1]+frecuency+'}}$'
              
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
        ax.set_title(name_statistics,fontsize=18)
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

def disaggregation_rainfall(x_series, y_series):
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


class Analysis(object):
    """"This function allows different analyses of the results obtained with respect to the calibration and actual data.
    Parameters
    ----------
    hiperparams : Object. 
        Object with calibration hyperparameters."""
    
    def __init__(self,CAL,SIM):
        self.CAL = CAL
        self.SIM = SIM
        self.figures = []
        self.names_figures = []
        
    def disaggregate_rainfall(self,x_series, y_series):
        self.hourly_disaggregation = disaggregation_rainfall(x_series,y_series)
        
    def compare_statistics_fig(self):
        self.figures.append(compare_statistics(self.CAL, self.SIM, self.SIM.temporal_resolution))
        self.names_figures.append('Statistical_comparison')
        
    def exceedence_probability_fig(self,Serie_Observed, Serie_Simulated):
        self.figures.append(exceedence_probability(Serie_Observed.mean(axis=1), Serie_Simulated.mean(axis=1), self.SIM.temporal_resolution))
        self.names_figures.append('Exceedence_probability')
        
    def figure_correlation_fig(self):
        [figures,names] = figure_correlation(self.CAL,self.SIM)
        
        for f in range(0, len(names)):
            self.figures.append(figures[f])
            self.names_figures.append(names[f]+'_analysis')
            
    def save_figures(self,path):
        for i,fig in enumerate(self.figures):
            fig.savefig(path+self.names_figures[i]+'.png',bbox_inches='tight')
        
        
    

