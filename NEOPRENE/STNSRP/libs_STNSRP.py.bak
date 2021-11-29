'''
Library containing the necessary functions to simulate the Neyman-Scott process. 
Rectangular Pulses (NSRP).

    Authors: 
        + Javier Díez Sierra
        + Salvador Navas Fernández
        + Manuel del Jesus
    
'''

'''Functions of the NSRP mode'''

from STNSRP.MathematicalPropertiesSTNSRP import *
from STNSRP.utils import *
import STNSRP.PSOMdJ as pso
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


def XI_MONTHS(Data, Dataframe_parametros, process, Seasonality):
    ##Calculo el parámetro de escala en cada localización
    Dataframe_xi=pd.DataFrame(index=Seasonality)
    for j in Data.columns:
        xis=list()
        for i in Seasonality: #Siempre ajusto el xi para los 12 meses
            posi=np.where(np.in1d(Data.index.month, i, assume_unique=True))[0]
            media_rain=np.nanmean(Data[j][posi])
            col=FIND_INDEX_MONTH(i, Dataframe_parametros)
            col=col[0]
            if process=='normal':
                if index_string(Dataframe_parametros.index, 'alpha')!=None: alpha=Dataframe_parametros[col]['alpha']
                else: alpha=1
                xis.append(scale_funtion(media_rain,Dataframe_parametros[col]['landa'],\
                                         Dataframe_parametros[col]['ipsilon'], Dataframe_parametros[col]['eta'],alpha, I_F))
            elif process=='cells':
                if index_string(Dataframe_parametros.index, 'alpha1')!=None: 
                    alpha1=Dataframe_parametros[col]['alpha1']; alpha1=Dataframe_parametros[col]['alpha2']; 
                else: alpha1=1; alpha2=1
                xi1=(scale_funtion(media_rain,Dataframe_parametros[col]['landa'],\
                         Dataframe_parametros[col]['ipsilon'], Dataframe_parametros[col]['eta1'],alpha1, I_F))
                xi2=(scale_funtion(media_rain,Dataframe_parametros[col]['landa'],\
                         Dataframe_parametros[col]['ipsilon'], Dataframe_parametros[col]['eta2'],alpha2, I_F))
                xis.append(xi1*media_rain,Dataframe_parametros[col]['alpha_p1']+\
                           xi2*(1-media_rain,Dataframe_parametros[col]['alpha_p1']))
            elif process=='storms':
                if index_string(Dataframe_parametros.index, 'alpha1')!=None: 
                    alpha1=Dataframe_parametros[col]['alpha1']; alpha1=Dataframe_parametros[col]['alpha2']; 
                else: alpha1=1; alpha2=1
                xi1=(scale_funtion(media_rain,Dataframe_parametros[col]['landa1'],\
                         Dataframe_parametros[col]['ipsilon1'], Dataframe_parametros[i]['eta1'],alpha1, I_F))
                xi2=(scale_funtion(media_rain,Dataframe_parametros[col]['landa1'],\
                         Dataframe_parametros[col]['ipsilon1'], Dataframe_parametros[i]['eta2'],alpha2, I_F))
                xis.append(xi1+xi2)

        Dataframe_xi[j]=xis

    ##Ordeno    
    Dataframe_xi_meses=pd.DataFrame(columns=Dataframe_xi.columns)
    for i in Dataframe_xi.index:
        if np.size(i)==1: Dataframe_xi_meses.loc[i]=Dataframe_xi.loc[[i]].values[0]
        else: 
            for mm in np.arange(1,13): 
                if mm in i: Dataframe_xi_meses.loc[mm]=Dataframe_xi.loc[[i]].values[0]

    #print(Dataframe_xi_meses)
    #Dataframe_xi_meses=Dataframe_xi_meses.sort()
    
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






