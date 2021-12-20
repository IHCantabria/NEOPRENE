'''
Library containing the necessary functions to simulate the Neyman-Scott process. 
Rectangular Pulses (NSRP).

    Authors: 
        + Javier Díez Sierra
        + Salvador Navas Fernández
        + Manuel del Jesus
    
'''

'''Functions of the NSRP mode'''

from NEOPRENE.STNSRP.MathematicalPropertiesSTNSRP import *
from NEOPRENE.STNSRP.utils import *
import NEOPRENE.STNSRP.PSOMdJ as pso
from haversine import haversine, Unit
import numpy as np
import pandas as pd
import scipy as sp
import math as mt
import datetime
import time
from scipy.interpolate import interp1d
from datetime import date
from datetime import timedelta
from scipy import stats
from scipy.optimize import curve_fit
import math as mt
from math import *
import tqdm

def scale_funtion(x, landa, ipsilon, eta, alpha):
    return (((landa*mt.gamma(2)*ipsilon))/(x*eta))

def index_string(param, string):
    if string in param:
        param_list=list(param)
        return param_list.index(string)
    else:
        return None

def calculate_statistics(Data,statistics,temporal_resolution):
    """ 
    Calculation of the statistics for a given station
    The statistics that can be calculated are:
    statistics=[mean, var_h, autocorr_l_h, fih_h, fiWW_h, fiDD_h, M3_h]
    (l=lag, h=aggregation levels)

    """

    statististics_values_real=list()
    if temporal_resolution=='d':
        t='D'
    elif temporal_resolution=='h':
        t='h'

    if np.sum(['var' in i for i in statistics])>=1:
        pos=np.where(['var' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_")[1])
            aux=Data.resample(str(h) + t).agg(pd.Series.sum, min_count=1); 
            statististics_values_real.append(np.sqrt(np.nanvar(aux))/np.nanmean(aux))
    if np.sum(['autocorr' in i for i in statistics])>=1:
        pos=np.where(['autocorr' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            l=int(statistics[ii].split("_")[1])
            h=int(statistics[ii].split("_")[2])
            aux=Data.resample(str(h) + t).agg(pd.Series.sum, min_count=1); 
            Autocorrelation_aux=aux.autocorr(lag=l) 
            if np.size(Autocorrelation_aux)>1: Autocorrelation_aux=Autocorrelation_aux[0]

            statististics_values_real.append(Autocorrelation_aux)
    if np.sum(['fih' in i for i in statistics])>=1:
        pos=np.where(['fih' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_")[1])
            statististics_values_real.append(fi_h(Data, h))
    if np.sum(['fiWW' in i for i in statistics])>=1:
        pos=np.where(['fiWW' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_")[1])
            statististics_values_real.append(fi_WW(Data, h))
    if np.sum(['fiDD' in i for i in statistics])>=1:
        pos=np.where(['fiDD' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_")[1])
            statististics_values_real.append(fi_DD(Data, h))
    if np.sum(['M3' in i for i in statistics])>=1:
        pos=np.where(['M3' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_")[1])
            aux=Data.resample(str(h) + t ).agg(pd.Series.sum, min_count=1);
            statististics_values_real.append((sp.stats.moment(aux, moment=3, nan_policy='omit'))/(np.nanvar(aux)**(3/2)))

    return statististics_values_real
    

def cross_corr_stationality_f(time_series, Seasonality, Attributes, func, coordinates, cross_corr_h, temporal_resolution):
    
    cross_corr_stationality={}
    
    cross_corr_dist_stationality={}
    

    for i in Seasonality:
        Data=time_series.copy()
        cross_corr=pd.DataFrame()

        m_T=list()
        for n in np.arange(1, 13): 
            if not(np.sum(i==n)>=1): m_T.append(n)
        if len(m_T)==0: Data_month=Data
        elif len(m_T)==1: 
            Data_month=Data
            Data_month[Data_month.index.month==m_T]=np.nan
        else:
            Data_month=Data; 
            for m_t in m_T:
                Data_month[Data_month.index.month==m_t]=np.nan
                
        h=int(cross_corr_h.split("_")[1])
        aux_month=Data.resample(str(h) + temporal_resolution).agg(pd.Series.sum, min_count=1); 

        cross_corr['dist'], cross_corr['cross_corr'] = cross_correlation(Attributes, aux_month, func, 11, coordinates)
        
        cross_corr_stationality[i] = cross_corr
        cross_corr_dist_stationality[i] = cross_correlation_dist(Attributes, aux_month, func, coordinates)

    return cross_corr_stationality, cross_corr_dist_stationality

def cross_correlation_dist(Stations, Series, funcion, coordinates):
    Spatial_correlation=pd.DataFrame(index=Series.columns, columns=Series.columns)
    Distance=pd.DataFrame(index=Series.columns, columns=Series.columns)

    for i, ii in enumerate(Spatial_correlation.index):
        x1=Stations[Stations['ID']==ii].X.values; x1=x1[0]
        y1=Stations[Stations['ID']==ii].Y.values; y1=y1[0]
        data1=Series[ii].values
        for j, jj, in enumerate(Spatial_correlation.columns):
            x2=Stations[Stations['ID']==jj].X.values; x2=x2[0]
            y2=Stations[Stations['ID']==jj].Y.values; y2=y2[0]
            if (coordinates.lower()=='geographical'):
                Distance[ii][jj]=haversine((y1, x1), (y2, x2), unit=Unit.KILOMETERS)
            elif coordinates.lower()=='utm':
                Distance[ii][jj]=distance_f(x1, y1, x2, y2)/1000
            else:
                raise Exception ('coordinates is wrongly defined, it should be geographical or utm')

                
            pos_no_nans=np.intersect1d(np.where(~np.isnan(Series[ii].values)), np.where(~np.isnan(Series[jj].values)))
            corr_s=stats.pearsonr(Series[ii].values[pos_no_nans], Series[jj].values[pos_no_nans])
            Spatial_correlation[ii][jj]=corr_s[0]

    Distance_correlation=pd.DataFrame()
    Distance_correlation['Corr']=np.reshape(Spatial_correlation.values, np.size(Spatial_correlation.values))
    Distance_correlation['dist']=np.reshape(Distance.values, np.size(Distance.values))
    
    Distance_correlation=pd.DataFrame(np.vstack({tuple(row) for row in Distance_correlation.values}), columns=['Corr', 'dist'])

    Distance_correlation = Distance_correlation[np.isfinite(Distance_correlation['Corr'])]
    
    return Distance_correlation


def cross_correlation(Stations, Series, funcion, divisions, coordinates):
    Spatial_correlation=pd.DataFrame(index=Series.columns, columns=Series.columns)
    Distance=pd.DataFrame(index=Series.columns, columns=Series.columns)

    for i, ii in enumerate(Spatial_correlation.index):
        x1=Stations[Stations['ID']==ii].X.values; x1=x1[0]
        y1=Stations[Stations['ID']==ii].Y.values; y1=y1[0]
        data1=Series[ii].values
        for j, jj, in enumerate(Spatial_correlation.columns):
            x2=Stations[Stations['ID']==jj].X.values; x2=x2[0]
            y2=Stations[Stations['ID']==jj].Y.values; y2=y2[0]
            if (coordinates.lower()=='geographical'):
                Distance[ii][jj]=haversine((y1, x1), (y2, x2), unit=Unit.KILOMETERS)
            elif coordinates.lower()=='utm':
                Distance[ii][jj]=distance_f(x1, y1, x2, y2)/1000
            else:
                raise Exception ('coordinates is wrongly defined, it should be geographical or utm')

                
            pos_no_nans=np.intersect1d(np.where(~np.isnan(Series[ii].values)), np.where(~np.isnan(Series[jj].values)))
            corr_s=stats.pearsonr(Series[ii].values[pos_no_nans], Series[jj].values[pos_no_nans])
            Spatial_correlation[ii][jj]=corr_s[0]

    Distance_correlation=pd.DataFrame()
    Distance_correlation['Corr']=np.reshape(Spatial_correlation.values, np.size(Spatial_correlation.values))
    Distance_correlation['dist']=np.reshape(Distance.values, np.size(Distance.values))
    
    Distance_correlation=pd.DataFrame(np.vstack({tuple(row) for row in Distance_correlation.values}), columns=['Corr', 'dist'])

    Distance_correlation = Distance_correlation[np.isfinite(Distance_correlation['Corr'])]

    sections=np.linspace(0, Distance_correlation['dist'].max(), divisions)
    average_points=(sections[1:] + sections[:-1]) / 2
    Mean_correlation=list()
    for i in range(len(average_points)):
        pos1=np.where(Distance_correlation['dist']>=sections[i]) 
        pos2=np.where(Distance_correlation['dist']<sections[i+1])
        pos=np.intersect1d(pos1, pos2)
        Mean_correlation.append(np.mean(Distance_correlation['Corr'][pos]))    
    
    
    xdata = Distance_correlation['dist'].values
    ydata = Distance_correlation['Corr'].values
    
    cross_correlation.Distance_correlation = Distance_correlation
    
    popt, pcov = curve_fit(func, xdata, ydata, maxfev=1000)
    
    
    average_points
    Mean_correlation=func(average_points, popt[0], popt[1], popt[2])
    
    return average_points, Mean_correlation


def XI_MONTHS(Data, Df_parameters, process):
    """Calculate the scale parameter for every gauge and season"""
    Df_xi=pd.DataFrame(index=Df_parameters.columns)
    for gauge in Data.columns:
        xis=list()
        for season in Df_parameters.columns:
            posi=np.where(np.in1d(Data.index.month, season, assume_unique=True))[0]
            mean_rain=np.nanmean(Data[gauge][posi])
            if process=='normal':
                if index_string(Df_parameters, 'alpha')!=None: alpha=Df_parameters[season]['alpha']
                else: alpha=1
                xis.append(scale_funtion(mean_rain,Df_parameters[season]['landa'],\
                                         Df_parameters[season]['ipsilon'], Df_parameters[season]['eta'],alpha))
            elif process=='storms':
                if index_string(Df_parameters.index, 'alpha1')!=None: 
                    alpha1=Df_parameters[season]['alpha1']; alpha1=Df_parameters[season]['alpha2']; 
                else: alpha1=1; alpha2=1
                xi1=(scale_funtion(mean_rain,Df_parameters[season]['landa1'],\
                         Df_parameters[season]['ipsilon1'], Df_parameters[season]['eta1'],alpha1))
                xi2=(scale_funtion(mean_rain,Df_parameters[season]['landa2'],\
                         Df_parameters[season]['ipsilon2'], Df_parameters[season]['eta2'],alpha2))
                xis.append(xi1+xi2)

        Df_xi[gauge]=xis

    # Sorting dataframe 
    Df_xi_meses=pd.DataFrame(columns=Df_xi.columns)
    for i in Df_xi.index:
        if np.size(i)==1: Df_xi_meses.loc[i]=Df_xi.loc[[i]].values[0]
        else: 
            for mm in np.arange(1,13): 
                if mm in i: Df_xi_meses.loc[mm]=Df_xi.loc[[i]].values[0]
    
    return Df_xi_meses


class evaluateInd_PSO(object):
    def __init__(self, vector, weights, process, statistics, storm_radius, cross_corr):
        self.v = vector
        self.w = weights
        self.p = process
        self.s = statistics
        self.storm_radius = storm_radius
        self.cross_corr = cross_corr

    def __call__(self, ind):
    
        landa=list(); mu_c=list(); eta=list(); xi=list(); betha=list(); alpha=list(); alpha_p=list();
        fi_may=list(); fi_may_s=list();
        if self.p=='normal' and self.storm_radius==False:
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(1); betha.append(ind[3]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[4]) 

        elif self.p=='normal' and self.storm_radius==True:
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(1); betha.append(ind[3]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[4]); fi_may_s.append(ind[5])

        elif self.p=='storms' and self.storm_radius==False:
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(1); betha.append(ind[3]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[4])
            landa.append(ind[5]); mu_c.append(ind[6]); eta.append(ind[7]); xi.append(1); betha.append(ind[8]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[9])

        elif self.p=='storms' and self.storm_radius==True:
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(1); betha.append(ind[3]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[4]); fi_may_s.append(ind[5])
            landa.append(ind[6]); mu_c.append(ind[7]); eta.append(ind[8]); xi.append(1); betha.append(ind[9]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[10]); fi_may_s.append(ind[11])
            

        d_e={}    
        ##Variance
        if np.sum(['var' in i for i in self.s])>=1:
            pos=np.where(['var' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",2)[1])
                a = sqrt(NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p))\
                    /NSRP_mean_ST(h, landa, mu_c, eta, xi, alpha, alpha_p)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
                
        ##Autocorrelation
        if np.sum(['autocorr' in i for i in self.s])>=1:
            pos=np.where(['autocorr' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                l=int(self.s[ii].split("_",3)[1])
                h=int(self.s[ii].split("_",3)[2])
                a=NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                    /NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a

        ##fi_h
        if np.sum(['fih' in i for i in self.s])>=1:
            pos=np.where(['fih' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",2)[1])
                a=NSRP_pdry(h, landa, mu_c, eta, betha, alpha_p)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a

        ##fi_WW
        if np.sum(['fiWW' in i for i in self.s])>=1:
            pos=np.where(['fiWW' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",2)[1])
                a = NSRP_fi_WW(h, landa, mu_c, eta, betha, alpha_p)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a

        ##fi_DD
        if np.sum(['fiDD' in i for i in self.s])>=1:
            pos=np.where(['fiDD' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",2)[1])
                a = NSRP_fi_DD(h, landa, mu_c, eta, betha, alpha_p)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a

        ##M3
        if np.sum(['M3' in i for i in self.s])>=1:
            pos=np.where(['M3' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",2)[1])
                a=NSRP_moments_order_3('Poisson',h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                /NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)**(3/2)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a

        ##Spatial correlation
        if self.storm_radius==False:
            if np.sum(['cross' in i for i in self.s])>=1:
                pos=np.where(['cross' in i for i in self.s]); pos=pos[0]
                for i, ii in enumerate(pos):
                    cross_corr_aux=self.cross_corr[(self.s[pos[i]])]['cross_corr']
                    dist_aux=self.cross_corr[(self.s[pos[i]])]['dist']
                    e_spatial_correlation=list()
                    for jj in range(len(cross_corr_aux)):
                        w1=self.w[ii]/11 ; l=0
                        h=int(self.s[ii].split("_")[1])
                        d=dist_aux[jj] 
                        a=NSRP_cross_correlation(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p, fi_may, d)\
                        /NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                        e_aux =w1*((1-(a/cross_corr_aux[jj]))**2 + (1-(cross_corr_aux[jj]/a))**2)
                        e_spatial_correlation.append(e_aux)
                        d_e['cross_corr_' + str(h) + '_' +  str(jj)]=e_aux

        if self.storm_radius==True:
            if np.sum(['cross' in i for i in self.s])>=1:
                pos=np.where(['cross' in i for i in self.s]); pos=pos[0]
                for i, ii in enumerate(pos):
                    cross_corr_aux=self.cross_corr[(self.s[pos[i]])]['cross_corr']
                    dist_aux=self.cross_corr[(self.s[pos[i]])]['dist']
                    e_spatial_correlation=list()
                    for jj in range(len(cross_corr_aux)):
                        w1=self.w[ii]/11; l=0
                        h=int(self.s[ii].split("_")[1])
                        d=dist_aux[jj] 
                        a=NSRP_cross_correlation_with_storm_radius(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p, \
                                                                       fi_may, d, fi_may_s)\
                        /NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                        e_aux =w1*((1-(a/cross_corr_aux[jj]))**2 + (1-(cross_corr_aux[jj]/a))**2)
                        e_spatial_correlation.append(e_aux)
                        d_e['cross_corr_' + str(h) + '_' +  str(jj)]=e_aux
        
        for k in d_e.keys():
            if isnan(d_e[k]): d_e[k] = 10000

        errors=list(); errors.append([d_e[i] for i in d_e.keys()])
        return np.sum(errors)

    def totalError(self, ind):

        landa=list(); mu_c=list(); eta=list(); xi=list(); betha=list(); alpha=list(); alpha_p=list();
        fi_may=list(); fi_may_s=list();
        if self.p=='normal' and self.storm_radius==False:
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(1); betha.append(ind[3]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[4]) 

        elif self.p=='normal' and self.storm_radius==True:
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(1); betha.append(ind[3]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[4]); fi_may_s.append(ind[5])

        elif self.p=='storms' and self.storm_radius==False:
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(1); betha.append(ind[3]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[4])
            landa.append(ind[5]); mu_c.append(ind[6]); eta.append(ind[7]); xi.append(1); betha.append(ind[8]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[9])

        elif self.p=='storms' and self.storm_radius==True:
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(1); betha.append(ind[3]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[4]); fi_may_s.append(ind[5])
            landa.append(ind[6]); mu_c.append(ind[7]); eta.append(ind[8]); xi.append(1); betha.append(ind[9]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[10]); fi_may_s.append(ind[11])
        
        d_e={}    
        ##Variance
        if np.sum(['var' in i for i in self.s])>=1:
            pos=np.where(['var' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = sqrt(NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p))\
                    /NSRP_mean_ST(h, landa, mu_c, eta, xi, alpha, alpha_p)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        ##Autocorrelation
        if np.sum(['autocorr' in i for i in self.s])>=1:
            pos=np.where(['autocorr' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                l=int(self.s[ii].split("_",3)[1])
                h=int(self.s[ii].split("_",3)[2])
                a=NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                    /NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        ##fi_h
        if np.sum(['fih' in i for i in self.s])>=1:
            pos=np.where(['fih' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a=NSRP_pdry(h, landa, mu_c, eta, betha, alpha_p)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        ##fi_WW
        if np.sum(['fiWW' in i for i in self.s])>=1:
            pos=np.where(['fiWW' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_fi_WW(h, landa, mu_c, eta, betha, alpha_p)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        ##fi_DD
        if np.sum(['fiDD' in i for i in self.s])>=1:
            pos=np.where(['fiDD' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_fi_DD(h, landa, mu_c, eta, betha, alpha_p)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        ##M3
        if np.sum(['M3' in i for i in self.s])>=1:
            pos=np.where(['M3' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a=NSRP_moments_order_3('Poisson',h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                /NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)**(3/2)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        ##Spatial correlation
        if self.storm_radius==False:
            if np.sum(['cross' in i for i in self.s])>=1:
                pos=np.where(['cross' in i for i in self.s]); pos=pos[0]
                for i, ii in enumerate(pos):
                    cross_corr_aux=self.cross_corr[(self.s[pos[i]])]['cross_corr']
                    dist_aux=self.cross_corr[(self.s[pos[i]])]['dist']
                    e_spatial_correlation=list()
                    for jj in range(len(cross_corr_aux)):
                        w1=self.w[ii]/11; l=0
                        #h=int(self.s[ii].split("_",3)[2])
                        h=int(self.s[ii].split("_")[1])
                        d=dist_aux[jj] 
                        a=NSRP_cross_correlation(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p, fi_may, d)\
                        /NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                        e_aux =w1*((1-(a/cross_corr_aux[jj]))**2 + (1-(cross_corr_aux[jj]/a))**2)
                        e_spatial_correlation.append(e_aux)
                        d_e['cross_corr_' + str(h) + '_' +  str(jj)]=e_aux
        if self.storm_radius==True:
            if np.sum(['cross' in i for i in self.s])>=1:
                pos=np.where(['cross' in i for i in self.s]); pos=pos[0]
                for i, ii in enumerate(pos):
                    cross_corr_aux=self.cross_corr[(self.s[pos[i]])]['cross_corr']
                    dist_aux=self.cross_corr[(self.s[pos[i]])]['dist']
                    e_spatial_correlation=list()
                    for jj in range(len(cross_corr_aux)):
                        w1=self.w[ii]/11; l=0
                        h=int(self.s[ii].split("_")[1])
                        d=dist_aux[jj] 
                        a=NSRP_cross_correlation_with_storm_radius(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p, \
                                                                       fi_may, d, fi_may_s)\
                        /NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                        e_aux =w1*((1-(a/cross_corr_aux[jj]))**2 + (1-(cross_corr_aux[jj]/a))**2)
                        e_spatial_correlation.append(e_aux)
                        d_e['cross_corr_' + str(h) + '_' +  str(jj)]=e_aux
        
        for k in d_e.keys():
            if isnan(d_e[k]): d_e[k] = 10000
                
        errors=list(); errors.append([d_e[i] for i in d_e.keys()])
        
        return np.sum(errors)
        
    def compute(self, ind):
        landa=list(); mu_c=list(); eta=list(); xi=list(); betha=list(); alpha=list(); alpha_p=list();
        fi_may=list(); fi_may_s=list();
        if self.p=='normal' and self.storm_radius==False:
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(1); betha.append(ind[3]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[4]) 

        elif self.p=='normal' and self.storm_radius==True:
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(1); betha.append(ind[3]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[4]); fi_may_s.append(ind[5])

        elif self.p=='storms' and self.storm_radius==False:
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(1); betha.append(ind[3]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[4])
            landa.append(ind[5]); mu_c.append(ind[6]); eta.append(ind[7]); xi.append(1); betha.append(ind[8]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[9])

        elif self.p=='storms' and self.storm_radius==True:
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(1); betha.append(ind[3]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[4]); fi_may_s.append(ind[5])
            landa.append(ind[6]); mu_c.append(ind[7]); eta.append(ind[8]); xi.append(1); betha.append(ind[9]);
            alpha.append(1); alpha_p.append(1); fi_may.append(ind[10]); fi_may_s.append(ind[11])
        
        v={}
        ##Variance
        if np.sum(['var' in i for i in self.s])>=1:
            pos=np.where(['var' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = sqrt(NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p))\
                    /NSRP_mean_ST(h, landa, mu_c, eta, xi, alpha, alpha_p)
                v['v' + str(ii)]=a
        ##Autocorrelation
        if np.sum(['autocorr' in i for i in self.s])>=1:
            pos=np.where(['autocorr' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                l=int(self.s[ii].split("_",3)[1])
                h=int(self.s[ii].split("_",3)[2])
                a=NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                    /NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                v['v' + str(ii)]=a
        ##fi_h
        if np.sum(['fih' in i for i in self.s])>=1:
            pos=np.where(['fih' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a=NSRP_pdry(h, landa, mu_c, eta, betha, alpha_p)
                v['v' + str(ii)]=a
        ##fi_WW
        if np.sum(['fiWW' in i for i in self.s])>=1:
            pos=np.where(['fiWW' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_fi_WW(h, landa, mu_c, eta, betha, alpha_p)
                v['v' + str(ii)]=a
        ##fi_DD
        if np.sum(['fiDD' in i for i in self.s])>=1:
            pos=np.where(['fiDD' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_fi_DD(h, landa, mu_c, eta, betha, alpha_p)
                v['v' + str(ii)]=a
        ##M3
        if np.sum(['M3' in i for i in self.s])>=1:
            pos=np.where(['M3' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a=NSRP_moments_order_3('Poisson',h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                /NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)**(3/2)
                v['v' + str(ii)]=a
        ##Spatial correlation
        if self.storm_radius==False:
            if np.sum(['cross' in i for i in self.s])>=1:
                pos=np.where(['cross' in i for i in self.s]); pos=pos[0]
                for i, ii in enumerate(pos):
                    cross_corr_aux=self.cross_corr[(self.s[pos[i]])]['cross_corr']
                    dist_aux=self.cross_corr[(self.s[pos[i]])]['dist']
                    e_spatial_correlation=list()
                    for jj in range(len(cross_corr_aux)):
                        w1=self.w[ii]/11; l=0
                        #h=int(self.s[ii].split("_",3)[2])
                        h=int(self.s[ii].split("_")[1])
                        d=dist_aux[jj] 
                        a=NSRP_cross_correlation(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p, fi_may, d)\
                        /NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                        v['cross_corr_' + str(h) + '_' +  str(jj)]=a
        if self.storm_radius==True:
            if np.sum(['cross' in i for i in self.s])>=1:
                pos=np.where(['cross' in i for i in self.s]); pos=pos[0]
                for i, ii in enumerate(pos):
                    cross_corr_aux=self.cross_corr[(self.s[pos[i]])]['cross_corr']
                    dist_aux=self.cross_corr[(self.s[pos[i]])]['dist']
                    e_spatial_correlation=list()
                    for jj in range(len(cross_corr_aux)):
                        w1=self.w[ii]/11; l=0
                        #h=int(self.s[ii].split("_",3)[2])
                        h=int(self.s[ii].split("_")[1])
                        d=dist_aux[jj] 
                        a=NSRP_cross_correlation_with_storm_radius(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p, \
                                                                       fi_may, d, fi_may_s)\
                        /NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                        v['cross_corr_' + str(h) + '_'+ str(jj)]=a
                            
        values=list(); values.append([v[i] for i in v.keys()])

        return values
    
from shapely.geometry import Point
from shapely.geometry import Polygon

def mean_area_intersect_rectangle_circle(x_lim, y_lim, radius, plott):
    ## We calculate the mean area intersecting a rectangle (fixed) and a circle whose center is moving along
    ## of the rectangle (it would be integrating between x and y. but I do it numerically)
    ##x_lim and y_lim would be the corners of the rectangle.
    
    min_x=np.min(x_lim)
    max_x=np.max(x_lim)
    min_y=np.min(y_lim)
    max_y=np.max(y_lim)
    
    rectangle= Polygon([(min_x, min_y), (min_x, max_y), (max_x, max_y), (max_x, min_y), (min_x, min_y)])

    xx=np.linspace(x_lim[0], x_lim[1], 100)
    yy=np.linspace(y_lim[0], y_lim[1], 100)
    [XX, YY]=np.meshgrid(xx, yy)
    XX_reshape=XX.reshape(np.size(XX), 1)
    YY_reshape=YY.reshape(np.size(YY), 1)
    
    mean_area=list()
    for i in range(len(XX_reshape)):
        circle=Point(XX_reshape[i],YY_reshape[i]).buffer(radius)
        difference=rectangle.intersection(circle)
        
        if plott=='True':
        
            fig = pyplot.figure(1, figsize=(10, 5), dpi=90)
            # 1: valid polygon
            ax = fig.add_subplot(121)
            patch = PolygonPatch(rectangle, facecolor='b', edgecolor='b', alpha=0.5, zorder=2)
            ax.add_patch(patch)
            patch = PolygonPatch(circle, facecolor='r', edgecolor='r', alpha=0.5, zorder=2)
            ax.add_patch(patch)
            xrange = [-1, 5]
            yrange = [-1, 5]
            ax.set_xlim(*xrange)
            ax.set_ylim(*yrange)
            
            # 2: intersection
            ax = fig.add_subplot(122)
            patch = PolygonPatch(difference, facecolor='g', edgecolor='g', alpha=0.5, zorder=2)
            ax.add_patch(patch)
            xrange = [-1, 5]
            yrange = [-1, 5]
            ax.set_xlim(*xrange)
            ax.set_ylim(*yrange)

        
        mean_area.append(difference.area)
    return (np.mean(mean_area))

def distance_f(x1, y1, x2, y2):
    dist=(((x1-x2)**2) + ((y1-y2)**2))**0.5
    return dist

def STNSRP_simulation(Params_month,Df_xi , XX, YY, year_ini, year_fin, temporal_resolution, process,coordinates,storm_radius, Seasonality, stations):
    
    landa = Params_month[Params_month.index=='landa'].values[0]
        
    Storm_origin = 1/landa ##Poisson process

    ipsilon = Params_month[Params_month.index=='ipsilon'].values[0]

    Number_cells_per_storm = ipsilon-1 ##Random poisson mean

    if process=='normal':
        eta = Params_month[Params_month.index=='eta'].values[0]

        Duration_cell =1 /eta ##exponencial random

    betha=Params_month[Params_month.index=='betha'].values[0]

    Dist_cell_origin=1/betha ##exponencial random

    fi_may=Params_month[Params_month.index=='fi_may'].values[0]

    Radio=1/fi_may

    if storm_radius==True:
        fi_may_s = Params_month[Params_month.index=='fi_may_s'].values[0]

        Radio_s  = 1/fi_may_s

    if coordinates.lower() == 'utm':#Lo paso a km
        Degrees_window=np.percentile(np.random.exponential(scale=1/np.max(fi_may), size=100000000),98)*1000 
    elif coordinates.lower() == 'geographical':

        Degrees_window=np.percentile(np.random.exponential(scale=1/np.max(fi_may), size=100000000),98)/111
    else:
        raise Exception ('coordinates is wrongly defined, it should be geographical or utm')

    P1=[np.min(XX)-Degrees_window, np.min(YY)-Degrees_window]; #print(P1)
    P2=[np.max(XX)+Degrees_window, np.min(YY)-Degrees_window]; #print(P2)
    P3=[np.max(XX)+Degrees_window, np.max(YY)+Degrees_window]; #print(P3)
    P4=[np.min(XX)-Degrees_window, np.max(YY)+Degrees_window]; #print(P4)
    
    xp=[P1[0], P2[0], P3[0], P4[0]]; yp=[P1[1], P2[1], P3[1], P4[1]];
    print('Simulation corners: \n xp:' + str(xp) + '\n yp:' + str(yp)) 


    if coordinates.lower() == 'utm':
        Distance_xx_square_km=distance_f(P1[0], P1[1], P2[0], P2[1])/1000;
        Distance_yy_square_km=distance_f(P1[0], P1[1], P4[0], P4[1])/1000;
    elif coordinates.lower() == 'geographical':
        Distance_xx_square_km=haversine((P1[1], P1[0]), (P2[1], P2[0]), unit=Unit.KILOMETERS); 
        Distance_yy_square_km=haversine((P1[1], P1[0]), (P4[1], P4[0]), unit=Unit.KILOMETERS);
    else:
        raise Exception ('coordinates is wrongly defined, it should be geographical or utm')

    Simulation_area=Distance_xx_square_km*Distance_yy_square_km; print("Simulation area (km²):" + str(Simulation_area))
    Simulation_area_degrees=abs((P1[0]-P2[0])*(P1[1]-P3[1]))

    fi_min=STNSRP_fi_min(Number_cells_per_storm, fi_may); print('Cells per storm per km²: \n' + str(fi_min))
    mu_c_area=fi_min*Simulation_area; print("Cells per storm in my simulation area: \n" + str(mu_c_area))

    Number_cells_per_storm=mu_c_area
    
    if storm_radius==True:
        ##Storm radius
        ## When entering the radius of the 
        ## storm radius as the cells that are left out I have to include them I say that if my storm radius covers for 
        ## example half of my simulation area then I will have twice as many storms. The problem is that I need to know 
        ## the average area that my storm intersects within my rectangle. If the storm were infinitesimal then the area
        ## measured would be the area of my storm and if the radius were larger than the major diagonal of the rectangle my mean area would be ##series my area of the rectangle.
        ## would be my area of the rectangle and then I wouldn't have to divide it or multiply it by anything.
        x_lim=[np.min(xp), np.max(xp)]; y_lim=[np.min(yp), np.max(yp)]
        Area_mean_storm=list()
        for rs in Radio_s:
            if coordinates.lower() == 'utm':
                rs_units = rs*1000
            elif coordinates.lower() == 'geographical':
                rs_units = rs/111
            else:
                raise Exception ('coordinates is wrongly defined, it should be geographical or utm')
            Area_mean_storm.append(mean_area_intersect_rectangle_circle(x_lim, y_lim, rs_units, 'False'))
            
        Area_mean_storm_units = Area_mean_storm

        
        Storm_origin_new = []
        for nst, st_m in enumerate(np.array(Area_mean_storm_units)):
            if st_m>Simulation_area_degrees:
                Storm_origin_new.append(Storm_origin[nst])
            else:
                Storm_origin_new.append(Storm_origin[nst]*(Simulation_area_degrees/Area_mean_storm_units[nst]))
        Storm_origin = Storm_origin_new
    
    time_d=pd.period_range(start=str((year_ini)),  end=str((year_fin)), freq='D')
    time_h=pd.period_range(start=str((year_ini)),  end=str((year_fin)), freq='h')

    Df_sim_join_day  = pd.DataFrame(index=time_d, columns = stations)
    Df_sim_join_hour = pd.DataFrame(index=time_h, columns = stations)

    for perio, monthss in enumerate(Seasonality):
        n_days=np.sum(np.in1d(Df_sim_join_day.index.month, monthss))
        position_day=np.in1d(Df_sim_join_day.index.month, monthss)
        position_hour=np.in1d(Df_sim_join_hour.index.month, monthss)

        Df_sim_day_aux  = pd.period_range('1800', periods=n_days, freq='D')
        Df_sim_hour_aux = pd.period_range('1800', periods=n_days, freq='h')
        Df_sim_min_aux  = pd.period_range('1800', periods=n_days, freq='min')

        time_star=datetime.datetime.strptime(str(Df_sim_day_aux[0]), '%Y-%m-%d')
        time_end=datetime.datetime.strptime(str(Df_sim_day_aux[-1]), '%Y-%m-%d')
        time_lapso=time_star

        time_storm_origins=list()

        n=0

        if np.size(monthss)==1:
            monthss=[monthss]
        else:
            monthss=[perio+1]

        while time_lapso < time_end:
            s = np.random.exponential(Storm_origin[monthss[0]-1], 1)

            if temporal_resolution=='d':
                time_lapso=time_lapso +  datetime.timedelta(days=s[0])
            elif temporal_resolution=='h':
                time_lapso=time_lapso +  datetime.timedelta(hours=s[0])

            time_storm_origins.append(time_lapso)
            n=n+1

        n_storms=len(time_storm_origins)

        if storm_radius==True:    
            ##Storm radius
            Rand_01_x=[np.random.uniform(0, 1) for i in range(n_storms)]
            x_storms=np.array(Rand_01_x)*(abs(P1[0]-P2[0])) + P1[0]; 
            Rand_01_y=[np.random.uniform(0, 1) for i in range(n_storms)]
            y_storms=np.array(Rand_01_y)*(abs(P1[1]-P4[1])) + P1[1]; 
            Storm_radius=list()
            for i in (time_storm_origins):
                Storm_radius.append(np.random.exponential(scale=Radio_s[i.month-1],size=1)[0])

        print('Number of storms ' +  str(n_storms) + ' for the months ' + str(Seasonality[perio]))
        Number_cell_per_storm=list()
        for i, ii in enumerate(time_storm_origins):
            Number_cell_per_storm.append(1 + (np.random.poisson(Number_cells_per_storm[monthss[0]-1], 1)))
        print('Number of rain cells per storm ' + str(np.mean(Number_cell_per_storm)))

        time0 = time_star
        time_ini_cells=list()
        time_fin_cells=list()
        Intensity_cells=list()
        Duration_hours_cells=list()
        radio_cells=list()
        x_cells=list(); y_cells=list()    
        
        for i in tqdm.tqdm(range(n_storms)):
            time1=time_storm_origins[i]

            Rand_01_x=[np.random.uniform(0, 1) for i in range(Number_cell_per_storm[i][0])]
            Rand_x=np.array(Rand_01_x)*(abs(P1[0]-P2[0])) + P1[0];
            Rand_01_y=[np.random.uniform(0, 1) for i in range(Number_cell_per_storm[i][0])]
            Rand_y=np.array(Rand_01_y)*(abs(P1[1]-P4[1])) + P1[1]; 

            Dist_hours_cell_sim=np.random.exponential(scale=Dist_cell_origin[monthss[0]-1], size=Number_cell_per_storm[i])
            Duration_hours_cell_sim=np.random.exponential(scale=Duration_cell[monthss[0]-1], size=Number_cell_per_storm[i])

            for j in range(Number_cell_per_storm[i][0]):

                pos_teta=(np.argmin(np.sqrt((XX-Rand_x[j])**2 + (YY-Rand_y[j])**2)))
                name_station=Df_xi.columns[pos_teta]

                Intensity_cell_sim=np.random.exponential(scale=1/Df_xi[name_station][monthss[0]],\
                                                              size=1)

                if temporal_resolution=='d':
                    time1_cell=time1 + datetime.timedelta(days=Dist_hours_cell_sim[j])   
                    time2_cell=time1_cell + datetime.timedelta(days=Duration_hours_cell_sim[j])

                elif temporal_resolution=='h':
                    time1_cell=time1 + datetime.timedelta(hours=Dist_hours_cell_sim[j])
                    time2_cell=time1_cell + datetime.timedelta(hours=Duration_hours_cell_sim[j])

                radio_cells.append(np.random.exponential(scale=Radio[monthss[0]-1], size=1))
                time_ini_cells.append(time1_cell)
                time_fin_cells.append(time2_cell)
                x_cells.append(Rand_x[j]); y_cells.append(Rand_y[j]);
                Intensity_cells.append(Intensity_cell_sim)
                Duration_hours_cells.append(Duration_hours_cell_sim[j])    
                
                
        ################################################################################
        time_ini_cells=np.array(time_ini_cells)
        time_fin_cells=np.array(time_fin_cells)
        Intensity_cells=np.array(Intensity_cells)
        x_cells=np.array(x_cells); y_cells=np.array(y_cells)
        radio_cells=np.array(radio_cells)
        Duration_hours_cells=np.array(Duration_hours_cells)
        #################################################################################
        ############################################################################        
        tt = pd.period_range(start=Df_sim_day_aux[0],end=Df_sim_day_aux[-1], freq='min')

        tt_ordinal=tt.astype(np.int32)*60*10**9

        ############################################################################
        Anhos=list()
        for i, ii in enumerate(time_fin_cells):
            aux=ii; year_aux=aux.year
            Anhos.append(year_aux)

        In_dates=np.where(np.array(Anhos)<year_fin)

        time_ini_cellss=time_ini_cells[In_dates[0]]

        time_fin_cellss=time_fin_cells[In_dates[0]]

        Intensity_cellss=Intensity_cells[In_dates[0]]

        Duration_hours_cellss=Duration_hours_cells[In_dates[0]]

        Intensity_cells=Intensity_cells[In_dates[0]]

        x_cells=x_cells[In_dates[0]]; y_cells=y_cells[In_dates[0]]; 

        radio_cells=radio_cells[In_dates[0]];

        for nunu, rr in enumerate(stations):
            x_aux=XX[nunu]; y_aux=YY[nunu]
            cells_wet=list()
            for ccc in range(len(Intensity_cells)):
                if coordinates.lower()=='geographical':
                    dist_cell_grid = (haversine((y_aux, x_aux), (y_cells[ccc], x_cells[ccc]), unit=Unit.KILOMETERS))
                    if radio_cells[ccc]  > dist_cell_grid:
                        cells_wet.append(ccc)
                elif coordinates.lower()=='utm':
                    dist_cell_grid=(distance_f(x_aux, y_aux, x_cells[ccc], y_cells[ccc]))
                    if radio_cells[ccc]>dist_cell_grid/1000:
                        cells_wet.append(ccc)
                else:
                    raise Exception ('coordinates is wrongly defined, it should be geographical or utm')



            time_ini_cells_aux=time_ini_cells[cells_wet]
            time_fin_cells_aux=time_fin_cells[cells_wet]
            Intensity_cells_aux=Intensity_cells[cells_wet]
            Duration_cells_aux=Duration_hours_cells[cells_wet]
            x_cells_aux=x_cells[cells_wet]; y_cells_aux=y_cells[cells_wet]
            radio_cells_aux=radio_cells[cells_wet]

            if storm_radius==True:
                time_ini_cells_aux_storm=list();
                time_fin_cells_aux_storm=list();
                Intensity_cells_aux_storm=list();
                Duration_cells_aux_storm=list()
                x_cells_aux_storm=list(); y_cells_aux_storm=list();
                radio_cells_aux_storm=list()
                for s, ss in enumerate(time_storm_origins):
                    if s+1<len(time_storm_origins):
                        posi_cells_tormenta=(time_ini_cells_aux>time_storm_origins[s]) & (time_ini_cells_aux<time_storm_origins[s+1])
                        posi_cells_tormenta=np.where(posi_cells_tormenta)[0]
                    else:
                        posi_cells_tormenta=(time_ini_cells_aux>time_storm_origins[s])
                        posi_cells_tormenta=np.where(posi_cells_tormenta)[0]

                    for c in range(len(posi_cells_tormenta)):
                        if coordinates.lower()=='geographical':
                            dist_cell_storm=(haversine((y_storms[s], x_storms[s]), \
                                                             (y_cells_aux[posi_cells_tormenta[c]],\
                                                             x_cells_aux[posi_cells_tormenta[c]]), unit=Unit.KILOMETERS))
                        elif coordinates.lower()=='utm':
                            dist_cell_storm=(distance_f(x_storms[s], y_storms[s], x_cells_aux[posi_cells_tormenta[c]], y_cells_aux[posi_cells_tormenta[c]]))/1000
                        else:
                            raise Exception ('coordinates is wrongly defined, it should be geographical or utm')

                        if Storm_radius[s]+radio_cells_aux[posi_cells_tormenta[c]]>dist_cell_storm:
                            time_ini_cells_aux_storm.append(time_ini_cells_aux[posi_cells_tormenta[c]]);
                            time_fin_cells_aux_storm.append(time_fin_cells_aux[posi_cells_tormenta[c]])
                            Intensity_cells_aux_storm.append(Intensity_cells_aux[posi_cells_tormenta[c]])
                            Duration_cells_aux_storm.append(Duration_cells_aux[posi_cells_tormenta[c]])
                            x_cells_aux_storm.append(x_cells_aux[posi_cells_tormenta[c]])
                            y_cells_aux_storm.append(y_cells_aux[posi_cells_tormenta[c]])
                            radio_cells_aux_storm.append(radio_cells_aux[posi_cells_tormenta[c]][0])

                time_ini_cells_aux=np.array(time_ini_cells_aux_storm)
                time_fin_cells_aux=np.array(time_fin_cells_aux_storm)
                Intensity_cells_aux=np.array(Intensity_cells_aux_storm)

            #xixi=Df_xi[Df_xi.index==monthss[0]].values[0][nunu]
            #Intensity_cells_aux2=np.random.exponential(scale=1/xixi, size=len(Intensity_cells_aux))
            #Intensity_cells_aux=Intensity_cells_aux2.copy()


            t_ini = np.hstack([time_ini_cells_aux, time_fin_cells_aux])

            if process=='cells':
                i_ini = np.hstack([Intensity_cells_aux.T, -Intensity_cells_aux.T])
                i_ini = i_ini[0]
            else: 
                i_ini = np.hstack([Intensity_cells_aux.T, -Intensity_cells_aux.T])
                i_ini = i_ini[0]

            orden = np.argsort(t_ini)
            t = t_ini[orden]

            i = np.cumsum(i_ini[orden])

            i[i<0] = 0

            rain=i.copy()

            t_ordinal= pd.PeriodIndex(t.astype(str),freq='N').astype(np.int32)

            rainfall = interp1d(t_ordinal, rain, kind="zero", bounds_error=False, fill_value=0.)
            

            rrain = rainfall(tt_ordinal)

            Date=pd.DataFrame(index=tt)

            if temporal_resolution=='d':
                Date['Rain']=rrain/(60*24)
            elif temporal_resolution=='h':
                Date['Rain']=rrain/(60)

            df = Date.copy()   

            del Date


            df2 = pd.DataFrame(df.groupby([df.index.year,df.index.month,df.index.day])[['Rain']].sum().values,
                               index = pd.period_range(start=df.index[0],end=df.index[-1],freq='D'))
            df3 = pd.DataFrame(df.groupby([df.index.year,df.index.month,df.index.day,df.index.hour])[['Rain']].sum().values,
                               index = pd.period_range(start=df.index[0],end=df.index[-1],freq='h'))


            Df_sim_join_day[rr].iloc[np.where(position_day)[0][0:len(df2)]]   = df2.values.flatten()


            Df_sim_join_hour[rr].iloc[np.where(position_hour)[0][0:len(df3)]] = df3.values.flatten()

            del rainfall, t_ordinal, rrain, i, t, orden, t_ini, i_ini
            
    return Df_sim_join_day, Df_sim_join_hour
    






