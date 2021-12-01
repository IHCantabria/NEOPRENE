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
import numpy as np
import pandas as pd
import scipy as sp
import math as mt
import datetime
from scipy.interpolate import interp1d
from datetime import date
from datetime import timedelta
from scipy import stats
from scipy.optimize import curve_fit
import math as mt
from math import *

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
            #Temporal correlatiin is igual to the covariance divided between the standard desviation.
            #I calculate directly the temporal correlation. (Actually I must fit the covariance) 
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

    for i in Seasonality:
        Data=time_series.copy()
        cross_corr=pd.DataFrame()

        ##Calculo los meses que no estan del 1 al 12. Pongo nans a los meses en los que no quiero calcular los
        ##estadísticos
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
                
        #h=int(cross_corr_h.split("_",3)[2])
        h=int(cross_corr_h.split("_")[1])
        aux_month=Data.resample(str(h) + temporal_resolution).agg(pd.Series.sum, min_count=1); 
        #aux_month=Datos.resample(str(h) + temporal_resolution , how='sum')
        #cross_corr_stationality_f.aa=aux_month

        cross_corr['dist'], cross_corr['cross_corr']=cross_correlation(Attributes, aux_month, func, 11, coordinates)
        
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
            if (coordinates=='Geographical') or (coordinates=='geographical'):
                Distancia[ii][jj]=haversine(x1, y1, x2, y2)
            elif coordinates=='UTM':
                Distancia[ii][jj]=distancia_f(x1, y1, x2, y2)
                Distancia=Distancia/1000#Lo paso a kilometros
                
            pos_no_nans=np.intersect1d(np.where(~np.isnan(Series[ii].values)), np.where(~np.isnan(Series[jj].values)))
            corr_s=stats.pearsonr(Series[ii].values[pos_no_nans], Series[jj].values[pos_no_nans])
            Correlacion_espacial[ii][jj]=corr_s[0]

    Correlacion_distancia=pd.DataFrame()
    Correlacion_distancia['Corr']=np.reshape(Correlacion_espacial.values, np.size(Correlacion_espacial.values))
    Correlacion_distancia['dist']=np.reshape(Distancia.values, np.size(Distancia.values))
    
    ##Quito los repetidos (solo me quedo con una parte de la diagonal)
    Correlacion_distancia=pd.DataFrame(np.vstack({tuple(row) for row in Correlacion_distancia.values}), columns=['Corr', 'dist'])
    ##Quito filas = nan
    Correlacion_distancia = Correlacion_distancia[np.isfinite(Correlacion_distancia['Corr'])]

    #fig= plt.figure(figsize=(20,8))
    #plt.plot(Correlacion_distancia['dist'], Correlacion_distancia['Corr'].values, '.r')
    #legnd = ['', 'Correlacion media']
    #plt.ylabel('correlation')
    #plt.xlabel('distance[km]')
    #plt.title('Spatial Correlation funtion')   

    #print(Correlacion_distancia) 
    
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


def XI_MONTHS(Data, Dataframe_parametros, process):
    """Calculate the scale parameter for every gauge and season"""
    Dataframe_xi=pd.DataFrame(index=Dataframe_parametros.columns)
    for gauge in Data.columns:
        xis=list()
        for season in Dataframe_parametros.columns:
            posi=np.where(np.in1d(Data.index.month, season, assume_unique=True))[0]
            mean_rain=np.nanmean(Data[gauge][posi])
            if process=='normal':
                if index_string(Dataframe_parametros, 'alpha')!=None: alpha=Dataframe_parametros[season]['alpha']
                else: alpha=1
                xis.append(scale_funtion(mean_rain,Dataframe_parametros[season]['landa'],\
                                         Dataframe_parametros[season]['ipsilon'], Dataframe_parametros[season]['eta'],alpha))
            elif process=='storms':
                if index_string(Dataframe_parametros.index, 'alpha1')!=None: 
                    alpha1=Dataframe_parametros[season]['alpha1']; alpha1=Dataframe_parametros[season]['alpha2']; 
                else: alpha1=1; alpha2=1
                xi1=(scale_funtion(mean_rain,Dataframe_parametros[season]['landa1'],\
                         Dataframe_parametros[season]['ipsilon1'], Dataframe_parametros[i]['eta1'],alpha1))
                xi2=(scale_funtion(mean_rain,Dataframe_parametros[season]['landa1'],\
                         Dataframe_parametros[season]['ipsilon1'], Dataframe_parametros[i]['eta2'],alpha2))
                xis.append(xi1+xi2)

        Dataframe_xi[gauge]=xis

    # Sorting dataframe 
    Dataframe_xi_meses=pd.DataFrame(columns=Dataframe_xi.columns)
    for i in Dataframe_xi.index:
        if np.size(i)==1: Dataframe_xi_meses.loc[i]=Dataframe_xi.loc[[i]].values[0]
        else: 
            for mm in np.arange(1,13): 
                if mm in i: Dataframe_xi_meses.loc[mm]=Dataframe_xi.loc[[i]].values[0]
    
    return Dataframe_xi_meses


class evaluateInd_PSO(object):
    def __init__(self, vector, weights, process, statistics, storm_radius, cross_corr):
        self.v = vector
        self.w = weights
        self.p = process
        self.s = statistics
        self.storm_radius = storm_radius
        self.cross_corr = cross_corr
        #self.I_F = 1
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
        ##Varianza
        if np.sum(['var' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['var' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",2)[1])
                a = sqrt(NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p))\
                    /NSRP_mean_ST(h, landa, mu_c, eta, xi, alpha, alpha_p)
                #a = self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
                
        ##Autocorrelacion
        if np.sum(['autocorr' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['autocorr' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                l=int(self.s[ii].split("_",3)[1])
                h=int(self.s[ii].split("_",3)[2])
                a=NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                    /NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                #a=self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a

        ##fi_h
        if np.sum(['fih' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['fih' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",2)[1])
                a=NSRP_pdry(h, landa, mu_c, eta, betha, alpha_p)
                #a = self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a

        ##fi_WW
        if np.sum(['fiWW' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['fiWW' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",2)[1])
                a = NSRP_fi_WW(h, landa, mu_c, eta, betha, alpha_p)
                #a = self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a

        ##fi_DD
        if np.sum(['fiDD' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['fiDD' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",2)[1])
                a = NSRP_fi_DD(h, landa, mu_c, eta, betha, alpha_p)
                #a = self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a

        ##M3
        if np.sum(['M3' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['M3' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",2)[1])
                a=NSRP_moments_order_3('Poisson',h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                /NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)**(3/2)
                #a=self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a

        ##Correlacion espacial
        if self.storm_radius==False:
            if np.sum(['cross' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
                pos=np.where(['cross' in i for i in self.s]); pos=pos[0]
                for i, ii in enumerate(pos):
                    cross_corr_aux=self.cross_corr[(self.s[pos[i]])]['cross_corr']
                    dist_aux=self.cross_corr[(self.s[pos[i]])]['dist']
                    e_correlacion_espacial=list()
                    for jj in range(len(cross_corr_aux)):
                        w1=1; l=0
                        #h=int(self.s[ii].split("_",3)[2])
                        h=int(self.s[ii].split("_")[1])
                        d=dist_aux[jj] 
                        a=NSRP_cross_correlation(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p, fi_may, d)\
                        /NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                        e_aux =w1*((1-(a/cross_corr_aux[jj]))**2 + (1-(cross_corr_aux[jj]/a))**2)
                        e_correlacion_espacial.append(e_aux)
                        d_e['cross_corr_' + str(h) + '_' +  str(jj)]=e_aux

        if self.storm_radius==True:
            if np.sum(['cross' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
                pos=np.where(['cross' in i for i in self.s]); pos=pos[0]
                for i, ii in enumerate(pos):
                    cross_corr_aux=self.cross_corr[(self.s[pos[i]])]['cross_corr']
                    dist_aux=self.cross_corr[(self.s[pos[i]])]['dist']
                    e_correlacion_espacial=list()
                    for jj in range(len(cross_corr_aux)):
                        w1=1; l=0
                        #h=int(self.s[ii].split("_",3)[2])
                        h=int(self.s[ii].split("_")[1])
                        d=dist_aux[jj] 
                        a=NSRP_cross_correlation_with_storm_radius(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p, \
                                                                       fi_may, d, fi_may_s)\
                        /NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                        e_aux =w1*((1-(a/cross_corr_aux[jj]))**2 + (1-(cross_corr_aux[jj]/a))**2)
                        e_correlacion_espacial.append(e_aux)
                        d_e['cross_corr_' + str(h) + '_' +  str(jj)]=e_aux
        
        for k in d_e.keys():
            if isnan(d_e[k]): d_e[k] = 10000

        errores=list(); errores.append([d_e[i] for i in d_e.keys()])
        
        return np.sum(errores)

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
        ##Varianza
        if np.sum(['var' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['var' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = sqrt(NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p))\
                    /NSRP_mean_ST(h, landa, mu_c, eta, xi, alpha, alpha_p)
                #a = self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        ##Autocorrelacion
        if np.sum(['autocorr' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['autocorr' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                l=int(self.s[ii].split("_",3)[1])
                h=int(self.s[ii].split("_",3)[2])
                a=NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                    /NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                #a = self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        ##fi_h
        if np.sum(['fih' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['fih' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a=NSRP_pdry(h, landa, mu_c, eta, betha, alpha_p)
                #a = self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        ##fi_WW
        if np.sum(['fiWW' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['fiWW' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_fi_WW(h, landa, mu_c, eta, betha, alpha_p)
                #a = self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        ##fi_DD
        if np.sum(['fiDD' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['fiDD' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_fi_DD(h, landa, mu_c, eta, betha, alpha_p)
                #a = self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        ##M3
        if np.sum(['M3' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['M3' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a=NSRP_moments_order_3('Poisson',h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                /NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)**(3/2)
                #a = self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        ##Correlacion espacial
        if self.storm_radius==False:
            if np.sum(['cross' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
                pos=np.where(['cross' in i for i in self.s]); pos=pos[0]
                for i, ii in enumerate(pos):
                    cross_corr_aux=self.cross_corr[(self.s[pos[i]])]['cross_corr']
                    dist_aux=self.cross_corr[(self.s[pos[i]])]['dist']
                    e_correlacion_espacial=list()
                    for jj in range(len(cross_corr_aux)):
                        w1=1; l=0
                        #h=int(self.s[ii].split("_",3)[2])
                        h=int(self.s[ii].split("_")[1])
                        d=dist_aux[jj] 
                        a=NSRP_cross_correlation(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p, fi_may, d)\
                        /NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                        e_aux =w1*((1-(a/cross_corr_aux[jj]))**2 + (1-(cross_corr_aux[jj]/a))**2)
                        e_correlacion_espacial.append(e_aux)
                        d_e['cross_corr_' + str(h) + '_' +  str(jj)]=e_aux
        if self.storm_radius==True:
            if np.sum(['cross' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
                pos=np.where(['cross' in i for i in self.s]); pos=pos[0]
                for i, ii in enumerate(pos):
                    cross_corr_aux=self.cross_corr[(self.s[pos[i]])]['cross_corr']
                    dist_aux=self.cross_corr[(self.s[pos[i]])]['dist']
                    e_correlacion_espacial=list()
                    for jj in range(len(cross_corr_aux)):
                        w1=1; l=0
                        #h=int(self.s[ii].split("_",3)[2])
                        h=int(self.s[ii].split("_")[1])
                        d=dist_aux[jj] 
                        a=NSRP_cross_correlation_with_storm_radius(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p, \
                                                                       fi_may, d, fi_may_s)\
                        /NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                        e_aux =w1*((1-(a/cross_corr_aux[jj]))**2 + (1-(cross_corr_aux[jj]/a))**2)
                        e_correlacion_espacial.append(e_aux)
                        d_e['cross_corr_' + str(h) + '_' +  str(jj)]=e_aux
        
        for k in d_e.keys():
            if isnan(d_e[k]): d_e[k] = 10000
                
        errores=list(); errores.append([d_e[i] for i in d_e.keys()])
        
        return np.sum(errores)
        
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
        ##Varianza
        if np.sum(['var' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['var' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = sqrt(NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p))\
                    /NSRP_mean_ST(h, landa, mu_c, eta, xi, alpha, alpha_p)
                v['v' + str(ii)]=a
        ##Autocorrelacion
        if np.sum(['autocorr' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['autocorr' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                l=int(self.s[ii].split("_",3)[1])
                h=int(self.s[ii].split("_",3)[2])
                a=NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                    /NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                v['v' + str(ii)]=a
        ##fi_h
        if np.sum(['fih' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['fih' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a=NSRP_pdry(h, landa, mu_c, eta, betha, alpha_p)
                v['v' + str(ii)]=a
        ##fi_WW
        if np.sum(['fiWW' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['fiWW' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_fi_WW(h, landa, mu_c, eta, betha, alpha_p)
                v['v' + str(ii)]=a
        ##fi_DD
        if np.sum(['fiDD' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['fiDD' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_fi_DD(h, landa, mu_c, eta, betha, alpha_p)
                v['v' + str(ii)]=a
        ##M3
        if np.sum(['M3' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['M3' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a=NSRP_moments_order_3('Poisson',h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                /NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)**(3/2)
                v['v' + str(ii)]=a
        ##Correlacion espacial
        if self.storm_radius==False:
            if np.sum(['cross' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
                pos=np.where(['cross' in i for i in self.s]); pos=pos[0]
                for i, ii in enumerate(pos):
                    cross_corr_aux=self.cross_corr[(self.s[pos[i]])]['cross_corr']
                    dist_aux=self.cross_corr[(self.s[pos[i]])]['dist']
                    e_correlacion_espacial=list()
                    for jj in range(len(cross_corr_aux)):
                        w1=1; l=0
                        #h=int(self.s[ii].split("_",3)[2])
                        h=int(self.s[ii].split("_")[1])
                        d=dist_aux[jj] 
                        a=NSRP_cross_correlation(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p, fi_may, d)\
                        /NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                        v['cross_corr_' + str(h) + '_' +  str(jj)]=a
        if self.storm_radius==True:
            if np.sum(['cross' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
                pos=np.where(['cross' in i for i in self.s]); pos=pos[0]
                for i, ii in enumerate(pos):
                    cross_corr_aux=self.cross_corr[(self.s[pos[i]])]['cross_corr']
                    dist_aux=self.cross_corr[(self.s[pos[i]])]['dist']
                    e_correlacion_espacial=list()
                    for jj in range(len(cross_corr_aux)):
                        w1=1; l=0
                        #h=int(self.s[ii].split("_",3)[2])
                        h=int(self.s[ii].split("_")[1])
                        d=dist_aux[jj] 
                        a=NSRP_cross_correlation_with_storm_radius(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p, \
                                                                       fi_may, d, fi_may_s)\
                        /NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                        v['cross_corr_' + str(h) + '_'+ str(jj)]=a
                            
        valores=list(); valores.append([v[i] for i in v.keys()])

        return valores
    
def NSRP_simulation(Params_month, year_ini, year_fin, temporal_resolution,process, Seasonality):

    landa=Params_month[Params_month.index=='landa'].values[0];Storm_origin=1/landa
    ipsilon=Params_month[Params_month.index=='ipsilon'].values[0];Number_cells_per_storm=ipsilon-1
    
    if process=='normal':
        eta=Params_month[Params_month.index=='eta'].values[0]; Time_cell=1/eta
        xi=Params_month[Params_month.index=='xi'].values[0]; Intensity_cell=1/xi

    betha=Params_month[Params_month.index=='betha'].values[0]; Dist_cell_origin=1/betha


    time_d=pd.period_range(start=str((year_ini)),  end=str((year_fin)), freq='D')
    time_h=pd.period_range(start=str((year_ini)),  end=str((year_fin)), freq='h')
    time_min=pd.period_range(start=str((year_ini)),  end=str((year_fin)), freq='min')

    Df_sim_join_day  =pd.DataFrame(index=time_d, columns = ['Rain'])
    Df_sim_join_hour =pd.DataFrame(index=time_h, columns = ['Rain'])


    Intensity_cells_total=list()
    Time_cells_total = list()

    for perio, monthss in enumerate(Seasonality):
        n_days=np.sum(np.in1d(Df_sim_join_day.index.month, monthss))
        position_day=np.in1d(Df_sim_join_day.index.month, monthss)
        position_hour=np.in1d(Df_sim_join_hour.index.month, monthss)

        Df_sim_day_aux  = pd.period_range('1800', periods=n_days, freq='D')
        Df_sim_hour_aux = pd.period_range('1800', periods=n_days*24, freq='h')
        Df_sim_min_aux  = pd.period_range('1800', periods=n_days*24*60, freq='min')

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
        
        Number_cell_per_storm=list()
        for i, ii in enumerate(time_storm_origins):
            Number_cell_per_storm.append(1 + (np.random.poisson(Number_cells_per_storm[monthss[0]-1], 1)))
        

        time0 = time_star
        time_ini_cells=list()
        time_fin_cells=list()
        Intensity_cells = list()
        Time_hours_cells = list()
        for i in range(n_storms):
            time1=time_storm_origins[i]
            Distance_hours_cell_sim=np.random.exponential(scale=Dist_cell_origin[monthss[0]-1], size=Number_cell_per_storm[i])

            if process=='normal':
                Time_hours_cell_sim=np.random.exponential(scale=Time_cell[monthss[0]-1], size=Number_cell_per_storm[i])
                
                Intensity_cell_sim=np.random.exponential(scale=Intensity_cell[monthss[0]-1], size=Number_cell_per_storm[i])

            for j in range(Number_cell_per_storm[i][0]):

                if temporal_resolution=='d':
                    time1_cell=time1 + datetime.timedelta(days=Distance_hours_cell_sim[j])
                    time2_cell=time1_cell + datetime.timedelta(days=Time_hours_cell_sim[j])

                elif temporal_resolution=='h':
                    time1_cell=time1 + datetime.timedelta(hours=Distance_hours_cell_sim[j])
                    time2_cell=time1_cell + datetime.timedelta(hours=Time_hours_cell_sim[j])

                Intensity=Intensity_cell_sim[j]
                time_ini_cells.append(time1_cell)
                time_fin_cells.append(time2_cell)
                Intensity_cells.append(Intensity)
                Time_hours_cells.append(Time_hours_cell_sim[j])

        ################################################################################
        time_ini_cells=np.array(time_ini_cells)
        time_fin_cells=np.array(time_fin_cells)
        Intensity_cells=np.array(Intensity_cells)
        Time_hours_cells=np.array(Time_hours_cells)
        #################################################################################
        ############################################################################
        tt=pd.period_range(start=Df_sim_day_aux[0],end=Df_sim_day_aux[-1], freq='min')
        tt_ordinal=tt.astype(np.int32)*60*10**9
        ############################################################################
        Anhos=list()
        for i, ii in enumerate(time_fin_cells):
            aux=ii; year_aux=aux.year
            Anhos.append(year_aux)
        In_datess=np.where(np.array(Anhos)<=Df_sim_day_aux.year[-1])
        time_ini_cellss=time_ini_cells[In_datess[0]]
        time_fin_cellss=time_fin_cells[In_datess[0]]
        Intensity_cellss=Intensity_cells[In_datess[0]]; Intensity_cells_total.append(Intensity_cellss)
        Time_hours_cellss=Time_hours_cells[In_datess[0]]; Time_cells_total.append(Time_hours_cellss)
        
        t_ini = np.hstack([time_ini_cellss, time_fin_cellss])
        if process=='cells':
            i_ini = np.hstack([Intensity_cellss.T, -Intensity_cellss.T])
            i_ini = i_ini[0]
        else: 
            i_ini = np.hstack([Intensity_cellss, -Intensity_cellss])
        orden = np.argsort(t_ini)
        t = t_ini[orden]
        i = np.cumsum(i_ini[orden])
        i[i<0] = 0
        rain=i.copy()
        t_ordinal = pd.PeriodIndex(t.astype(str),freq='N').astype(np.int32)
        if np.size(t_ordinal)==0:
            a=1
        else:
            rainfall = interp1d(t_ordinal, rain, kind="zero", bounds_error=False, fill_value=0.)
            rr = rainfall(tt_ordinal)
            Date=pd.DataFrame(index=tt)
            if temporal_resolution=='d':
                Date['Rain']=rr/(60*24)
            elif temporal_resolution=='h':
                Date['Rain']=rr/(60)
            df = Date.copy()
            del Date
            df2 = pd.DataFrame(df.groupby([df.index.year,df.index.month,df.index.day])[['Rain']].sum().values,index = pd.period_range(start=df.index[0],end=df.index[-1],freq='D'))
            df3 = pd.DataFrame(df.groupby([df.index.year,df.index.month,df.index.day,df.index.hour])[['Rain']].sum().values,index = pd.period_range(start=df.index[0],end=df.index[-1],freq='h'))


            Df_sim_join_day['Rain'].iloc[np.where(position_day)[0][0:len(df2)]]=df2.values.flatten()
            
            Df_sim_join_hour['Rain'].iloc[np.where(position_hour)[0][0:len(df3)]]=df3.values.flatten()
            
        del rr, rainfall, tt, t_ordinal,tt_ordinal, rain, i, t, orden, t_ini, i_ini,df, df3, df2
    return Df_sim_join_hour, Df_sim_join_day,\
           np.hstack(Intensity_cells_total), np.hstack(Time_cells_total)


from shapely.geometry import Point
from shapely.geometry import Polygon
#from descartes.patch import PolygonPatch
def mean_area_intersect_rectangle_circle(x_lim, y_lim, radius, plott):
    ##Calculo el aréa media que intersecta un rectangulo (fijo) y un circulo cuyo centro se va moviendo a lo largo
    ## del rectangulo (Sería integrar entre x e y. pero lo hago numericamente)
    ##x_lim y y_lim serían las esquinas del rectangulo
    
    min_x=np.min(x_lim)
    max_x=np.max(x_lim)
    min_y=np.min(y_lim)
    max_y=np.max(y_lim)
    
    rectangle= Polygon([(min_x, min_y), (min_x, max_y), (max_x, max_y), (max_x, min_y), (min_x, min_y)])##creo rectangulo

    xx=np.linspace(x_lim[0], x_lim[1], 100)
    yy=np.linspace(y_lim[0], y_lim[1], 100)
    [XX, YY]=np.meshgrid(xx, yy)
    XX_reshape=XX.reshape(np.size(XX), 1)
    YY_reshape=YY.reshape(np.size(YY), 1)
    
    area_medio=list()
    for i in range(len(XX_reshape)):
        circle=Point(XX_reshape[i],YY_reshape[i]).buffer(radius)##creo circulo
        diferencia=rectangle.intersection(circle)##calculo intersección
        
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
            #ax.set_xticks(range(*xrange) + [xrange[-1]])
            ax.set_ylim(*yrange)
            #ax.set_yticks(range(*yrange) + [yrange[-1]])
            #ax.set_aspect(1)

            # 2: intersection
            ax = fig.add_subplot(122)
            patch = PolygonPatch(diferencia, facecolor='g', edgecolor='g', alpha=0.5, zorder=2)
            ax.add_patch(patch)
            xrange = [-1, 5]
            yrange = [-1, 5]
            ax.set_xlim(*xrange)
            #ax.set_xticks(range(*xrange) + [xrange[-1]])
            ax.set_ylim(*yrange)
            #ax.set_yticks(range(*yrange) + [yrange[-1]])
            #ax.set_aspect(1)
        
        area_medio.append(diferencia.area)
    return (np.mean(area_medio))

def distancia_f(x1, y1, x2, y2):
    dist=((x1-x2)**2 + (y1-y2)**2)**0.5
    return dist

def STNSRP_simulation(Params_month, XX, YY, year_ini, year_fin, temporal_resolution, process,coordinates,storm_radius, Seasonality):
    
    
    
    landa = Params_month[Params_month.index=='landa'].values[0]
        
    Storm_origin = 1/landa##Poisson process

    ipsilon = Params_month[Params_month.index=='ipsilon'].values[0]

    Number_cells_per_storm = ipsilon-1##Random poisson mean ¿ipsilon-1?

    if process=='normal':
        eta = Params_month[Params_month.index=='eta'].values[0]

        Duracion_cell =1 /eta ##exponencial random

    betha=Params_month[Params_month.index=='betha'].values[0]

    Dist_cell_origin=1/betha ##exponencial random

    fi_may=Params_month[Params_month.index=='fi_may'].values[0]

    Radio=1/fi_may

    if storm_radius==True:
        fi_may_s = Params_month[Params_month.index=='fi_may_s'].values[0]

        Radio_s  = 1/fi_may_s

    if coordinates == 'UTM':#Lo paso a km

        Grados_ventana=np.max(np.random.exponential(scale=1/np.max(fi_may), size=100000000))*1000 #Lo paso a metros
    else:
        Grados_ventana=np.max(np.random.exponential(scale=1/np.max(fi_may), size=100000000))/111 #Lo paso a grados

    P1=[np.min(XX)-Grados_ventana, np.min(YY)-Grados_ventana]; print(P1)
    P2=[np.max(XX)+Grados_ventana, np.min(YY)-Grados_ventana]; print(P2)
    P3=[np.max(XX)+Grados_ventana, np.max(YY)+Grados_ventana]; print(P3)
    P4=[np.min(XX)-Grados_ventana, np.max(YY)+Grados_ventana]; print(P4)
    xp=[P1[0], P2[0], P3[0], P4[0]]; yp=[P1[1], P2[1], P3[1], P4[1]];


    if coordinates=='UTM':
        Distnacia_xx_cuadrado_km=distancia_f(P1[0], P1[1], P2[0], P2[1])/1000; print(Distnacia_xx_cuadrado_km)
        Distnacia_yy_cuadrado_km=distancia_f(P1[0], P1[1], P4[0], P4[1])/1000; print(Distnacia_yy_cuadrado_km)
    else:
        Distnacia_xx_cuadrado_km=haversine((P1[0], P1[1]), (P2[0], P2[1])); print(Distnacia_xx_cuadrado_km)
        Distnacia_yy_cuadrado_km=haversine((P1[0], P1[1]), (P4[0], P4[1])); print(Distnacia_yy_cuadrado_km)

    Area_simulacion=Distnacia_xx_cuadrado_km*Distnacia_yy_cuadrado_km; print(Area_simulacion)
    Area_simulacion_degrees=abs((P1[0]-P2[0])*(P1[1]-P3[1]))

    fi_min=STNSRP_fi_min(Number_cells_per_storm, fi_may); print(str(fi_min) +' Celdas por km² y por tormenta')#
    mu_c_area=fi_min*Area_simulacion; print(str(mu_c_area) + ' Celdas por tormenta en mi area de simulacion')

    Number_cells_per_storm=mu_c_area
    
    print('Storm ini = ' + str(Storm_origin))

    if storm_radius==True:
        ##Storm radius
        ## Al introducir el radio de la 
        ##tormenta como las celdas que quedan fuera tengo que incluirlas digo que si mi radio de tormenta cubre por
        ##ejemplo la mitad de mi area de simulación entonces tendre el doble de tormentas. El problema es que necesito conocer 
        ##el area media que ocupa intersecta mi tormenta dentro de mi rectangulo. Si la tormenta fuera infinitesimal el area
        ##medí sería el area de mi tormenta y si el radio fuera mas grande que la diagonal mayor del rectangulo mi area media
        ## serie mi area del rectangulo y entonces no habría que dividirlo ni multiplicarlo por nada.
        x_lim=[np.min(xp), np.max(xp)]; y_lim=[np.min(yp), np.max(yp)]
        Area_media_tormenta=list()
        for rs in Radio_s:
            Area_media_tormenta.append(mean_area_intersect_rectangle_circle(x_lim, y_lim, rs*1000, 'False'))
        print('Area simulacion degrees ' + str(Area_simulacion_degrees))
        print('Area media tormenta ' + str(Area_media_tormenta))
        Storm_origin_with_storm_radious=(np.array(Area_media_tormenta)/np.array(Area_simulacion_degrees))\
                                        *np.array(Storm_origin)
        Storm_origin=Storm_origin_with_storm_radious
    
    
    time_d=pd.period_range(start=str((year_ini)),  end=str((year_fin)), freq='D')
    time_h=pd.period_range(start=str((year_ini)),  end=str((year_fin)), freq='h')
    time_min=pd.period_range(start=str((year_ini)),  end=str((year_fin)), freq='min')

    Df_sim_join_day  = pd.DataFrame(index=time_d, columns = Datos_.columns)
    Df_sim_join_hour = pd.DataFrame(index=time_h, columns = Datos_.columns)
    Df_sim_join_min  = pd.DataFrame(index=time_min, columns = Datos_.columns)

    Intensidad_cells_total=list()
    Duracion_cells_total=list()

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
                time_lapso=time_lapso +  datetime.timedelta(days=s[0])##CAMBIAR dependeindo si estas en h o d
            elif temporal_resolution=='h':
                time_lapso=time_lapso +  datetime.timedelta(hours=s[0])##CAMBIAR dependeindo si estas en h o d

            time_storm_origins.append(time_lapso)
            n=n+1

        n_storms=len(time_storm_origins)

        if storm_radius==True:    
            ##Storm radius (Simulo la posición los centros de las tormentas)
            Rand_01_x=[np.random.uniform(0, 1) for i in range(n_storms)]
            x_storms=np.array(Rand_01_x)*(abs(P1[0]-P2[0])) + P1[0]; 
            Rand_01_y=[np.random.uniform(0, 1) for i in range(n_storms)]
            y_storms=np.array(Rand_01_y)*(abs(P1[1]-P4[1])) + P1[1]; 
            Storm_radius=list()
            for i in (time_storm_origins):
                Storm_radius.append(np.random.exponential(scale=Radio_s[i.month-1],size=1)[0])

        print('Numero de tormentas ' +  str(n_storms) + ' para los meses ' + str(monthss))
        Number_cell_per_storm=list()
        for i, ii in enumerate(time_storm_origins):
            Number_cell_per_storm.append(1 + (np.random.poisson(Number_cells_per_storm[monthss[0]-1], 1)))
        print('Numero de celdas de lluvia por tormenta ' + str(np.mean(Number_cell_per_storm)))

        time0 = time_star
        time_ini_cells=list()
        time_fin_cells=list()
        Intensidad_cells=list()
        Duracion_horas_cells=list()
        radio_cells=list()
        x_cells=list(); y_cells=list()    

        for i in range(n_storms):##numero de rain cells
            time1=time_storm_origins[i] #ojo horas!

            Rand_01_x=[np.random.uniform(0, 1) for i in range(Number_cell_per_storm[i][0])]
            Rand_x=np.array(Rand_01_x)*(abs(P1[0]-P2[0])) + P1[0];
            Rand_01_y=[np.random.uniform(0, 1) for i in range(Number_cell_per_storm[i][0])]
            Rand_y=np.array(Rand_01_y)*(abs(P1[1]-P4[1])) + P1[1]; 

            Distancia_horas_cell_sim=np.random.exponential(scale=Dist_cell_origin[monthss[0]-1], size=Number_cell_per_storm[i])
            Duracion_horas_cell_sim=np.random.exponential(scale=Duracion_cell[monthss[0]-1], size=Number_cell_per_storm[i])

            for j in range(Number_cell_per_storm[i][0]):#Nuevo

                if spatially_varying_intensity=='KNN':
                    pos_teta=(np.argmin(np.sqrt((Estaciones.X-Rand_x[j])**2 + (Estaciones.Y-Rand_y[j])**2)))
                    name_estacion=Dataframe_xi.columns[pos_teta]

                    Intensidad_cell_sim=np.random.exponential(scale=1/Dataframe_xi[name_estacion][monthss[0]],\
                                                                  size=1)

                if temporal_resolution=='d':
                    time1_cell=time1 + datetime.timedelta(days=Distancia_horas_cell_sim[j])##CAMBIAR dependeindo si estas en h o d     
                    time2_cell=time1_cell + datetime.timedelta(days=Duracion_horas_cell_sim[j])##CAMBIAR dependeindo si estas en h o d

                elif temporal_resolution=='h':
                    time1_cell=time1 + datetime.timedelta(hours=Distancia_horas_cell_sim[j])##CAMBIAR dependeindo si estas en h o d
                    time2_cell=time1_cell + datetime.timedelta(hours=Duracion_horas_cell_sim[j])##CAMBIAR dependeindo si estas en h o d

                #Intensidad.append(Intensidad_cell_sim)
                radio_cells.append(np.random.exponential(scale=Radio[monthss[0]-1], size=1))
                time_ini_cells.append(time1_cell)
                time_fin_cells.append(time2_cell)
                x_cells.append(Rand_x[j]); y_cells.append(Rand_y[j]);
                Intensidad_cells.append(Intensidad_cell_sim)
                Duracion_horas_cells.append(Duracion_horas_cell_sim[j])       
        ################################################################################
        time_ini_cells=np.array(time_ini_cells)
        time_fin_cells=np.array(time_fin_cells)
        Intensidad_cells=np.array(Intensidad_cells)
        x_cells=np.array(x_cells); y_cells=np.array(y_cells)
        radio_cells=np.array(radio_cells)
        Duracion_horas_cells=np.array(Duracion_horas_cells)
        #################################################################################
        ############################################################################        
        tt = pd.period_range(start=Df_sim_day_aux[0],end=Df_sim_day_aux[-1], freq='min')

        tt_ordinal=tt.astype(np.int32)*60*10**9

        ############################################################################
        Anhos=list()
        for i, ii in enumerate(time_fin_cells):
            aux=ii; year_aux=aux.year
            Anhos.append(year_aux)
        ##Quito las cendas que dentro de una tormenta caen fuera del año límite

        Dentro_fechass=np.where(np.array(Anhos)<year_fin)

        time_ini_cellss=time_ini_cells[Dentro_fechass[0]]

        time_fin_cellss=time_fin_cells[Dentro_fechass[0]]

        Intensidad_cellss=Intensidad_cells[Dentro_fechass[0]]

        Duracion_horas_cellss=Duracion_horas_cells[Dentro_fechass[0]]

        Intensidad_cells=Intensidad_cells[Dentro_fechass[0]]

        x_cells=x_cells[Dentro_fechass[0]]; y_cells=y_cells[Dentro_fechass[0]]; 

        radio_cells=radio_cells[Dentro_fechass[0]];

        for nunu, rr in enumerate(Datos_.columns):
            #time.sleep(0.01)
            ##Veo las celdas de lluvia que tocan el primer grid
            x_aux=XX[nunu]; y_aux=YY[nunu]
            celdas_mojan=list()
            for ccc in range(len(Intensidad_cells)):
                if coordinates=='geographical':
                    #distancia_celda_grid= (haversine(x_aux, y_aux, x_cells[ccc], y_cells[ccc]))
                    distancia_celda_grid = (haversine((x_aux, y_aux), (x_cells[ccc], y_cells[ccc])))
                    if radio_cells[ccc]  > distancia_celda_grid:
                        celdas_mojan.append(ccc)
                else:
                    distancia_celda_grid=(distancia_f(x_aux, y_aux, x_cells[ccc], y_cells[ccc]))
                    if radio_cells[ccc]>distancia_celda_grid/1000:
                        celdas_mojan.append(ccc)

            #zeros=np.zeros((len(Df_sim_join_min.index), 1))
            #aux_t=np.zeros((len(Df_sim_join_min.index), 1))

            time_ini_cells_aux=time_ini_cells[celdas_mojan]
            time_fin_cells_aux=time_fin_cells[celdas_mojan]
            Intensidad_cells_aux=Intensidad_cells[celdas_mojan]
            Duracion_cells_aux=Duracion_horas_cells[celdas_mojan]
            x_cells_aux=x_cells[celdas_mojan]; y_cells_aux=y_cells[celdas_mojan]
            radio_cells_aux=radio_cells[celdas_mojan]

            if storm_radius==True:
                ##Ahora voy a quitar las celdas que no mojan por que se salen del radio de la tormenta
                time_ini_cells_aux_storm=list();
                time_fin_cells_aux_storm=list();
                Intensidad_cells_aux_storm=list();
                Duracion_cells_aux_storm=list()
                x_cells_aux_storm=list(); y_cells_aux_storm=list();
                radio_cells_aux_storm=list()
                for s, ss in enumerate(time_storm_origins):
                    if s+1<len(time_storm_origins):
                        posi_celdas_tormenta=(time_ini_cells_aux>time_storm_origins[s]) & (time_ini_cells_aux<time_storm_origins[s+1])
                        posi_celdas_tormenta=np.where(posi_celdas_tormenta)[0]
                    else:
                        posi_celdas_tormenta=(time_ini_cells_aux>time_storm_origins[s])
                        posi_celdas_tormenta=np.where(posi_celdas_tormenta)[0]

                    for c in range(len(posi_celdas_tormenta)):
                        ##Veo si el pixel está dentro de la tormenta 
                        distancia_celda_storm=(haversine((x_storms[s], y_storms[s]), \
                                                         (x_cells_aux[posi_celdas_tormenta[c]],\
                                                         y_cells_aux[posi_celdas_tormenta[c]])))
                        if Storm_radius[s]+radio_cells_aux[posi_celdas_tormenta[c]]>distancia_celda_storm:
                            time_ini_cells_aux_storm.append(time_ini_cells_aux[posi_celdas_tormenta[c]]);
                            time_fin_cells_aux_storm.append(time_fin_cells_aux[posi_celdas_tormenta[c]])
                            Intensidad_cells_aux_storm.append(Intensidad_cells_aux[posi_celdas_tormenta[c]])
                            Duracion_cells_aux_storm.append(Duracion_cells_aux[posi_celdas_tormenta[c]])
                            x_cells_aux_storm.append(x_cells_aux[posi_celdas_tormenta[c]])
                            y_cells_aux_storm.append(y_cells_aux[posi_celdas_tormenta[c]])
                            radio_cells_aux_storm.append(radio_cells_aux[posi_celdas_tormenta[c]][0])

                time_ini_cells_aux=np.array(time_ini_cells_aux_storm)
                time_fin_cells_aux=np.array(time_fin_cells_aux_storm)
                Intensidad_cells_aux=np.array(Intensidad_cells_aux_storm)

            #if multisite==True:
            ###################OJO estoy haciendo simulación multisite. Cuando hago IDF o kNN cambia mucho la
            #precipitación media
            xixi=Dataframe_xi[Dataframe_xi.index==monthss[0]].values[0][nunu]
            Intensidad_cells_aux2=np.random.exponential(scale=1/xixi, size=len(Intensidad_cells_aux))
            Intensidad_cells_aux=Intensidad_cells_aux2.copy()
            ####################
            ##Acumulo cada celda de lluvia en minutos para luego agruparlo en horas

            t_ini = np.hstack([time_ini_cells_aux, time_fin_cells_aux])

            if process=='cells':
                i_ini = np.hstack([Intensidad_cells_aux.T, -Intensidad_cells_aux.T])
                i_ini = i_ini[0]
            else: 
                i_ini = np.hstack([Intensidad_cells_aux.T, -Intensidad_cells_aux.T])
                #if multisite==True:
                a=1
                #else:
                #    i_ini = i_ini[0]#OJO PONER

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
                Date['Rain']=rrain/(60)###calculo todo en horas y lo estoy interpolando a minutos por eso lo divido entre 60. Si interpolase

            df = Date.copy()   

            del Date


            df2 = pd.DataFrame(df.groupby([df.index.year,df.index.month,df.index.day])[['Rain']].sum().values,
                               index = pd.period_range(start=df.index[0],end=df.index[-1],freq='D'))
            df3 = pd.DataFrame(df.groupby([df.index.year,df.index.month,df.index.day,df.index.hour])[['Rain']].sum().values,
                               index = pd.period_range(start=df.index[0],end=df.index[-1],freq='h'))

            #Df_sim_join_day['Rain'].iloc[np.where(position_day)[0][0:len(df2)]]=df2.values.flatten()

            Df_sim_join_day[rr].iloc[np.where(position_day)[0][0:len(df2)]]   = df2.values.flatten()


            Df_sim_join_hour[rr].iloc[np.where(position_hour)[0][0:len(df3)]] = df3.values.flatten()

            del rainfall, t_ordinal, rrain, i, t, orden, t_ini, i_ini
            
    return Df_sim_join_day, Df_sim_join_hour
    






