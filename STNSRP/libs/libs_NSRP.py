'''
Librería que contine las funciones necesarias para simular el proceso Neyman-Scott 
Rectangular Pulses (NSRP).

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
from scipy.interpolate import interp1d
from datetime import date
from datetime import timedelta
from scipy import stats

def calculate_statistics(Datos,statistics,temporal_resolution):
    """Calculo los estaísticos para una estacion determinada
    Los estadísticos que se pueden calcular son:
                statistics=[mean, var_h, autocorr_l_h, fih_h, fiWW_h, fiDD_h, M3_h]
                        (l=lag, h=aggregation levels)

    """

    statististics_values_real=list()
    if temporal_resolution=='d':
        t='D'
    elif temporal_resolution=='h':
        t='h'

    if np.sum(['mean' in i for i in statistics])>=1:#Quiere decir que hay que ajustar ese estadístico
        statististics_values_real.append(np.nanmean(Datos))
    if np.sum(['var' in i for i in statistics])>=1:
        pos=np.where(['var' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",1)[1])
            aux=Datos.resample(str(h) + t , how='sum'); 
            statististics_values_real.append(np.nanvar(aux))
    if np.sum(['autocorr' in i for i in statistics])>=1:
        pos=np.where(['autocorr' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            l=int(statistics[ii].split("_",3)[1])
            h=int(statistics[ii].split("_",3)[2])
            aux=Datos.resample(str(h) + t , how='sum'); 
            Autocorrelation_aux=aux.autocorr(lag=l) 
            if np.size(Autocorrelation_aux)>1: Autocorrelation_aux=Autocorrelation_aux[0] 
            statististics_values_real.append(Autocorrelation_aux)
    if np.sum(['fih' in i for i in statistics])>=1:
        pos=np.where(['fih' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",1)[1])
            statististics_values_real.append(fi_h(Datos, h))
    if np.sum(['fiWW' in i for i in statistics])>=1:
        pos=np.where(['fiWW' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",1)[1])
            statististics_values_real.append(fi_WW(Datos, h))
    if np.sum(['fiDD' in i for i in statistics])>=1:
        pos=np.where(['fiDD' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",1)[1])
            statististics_values_real.append(fi_DD(Datos, h))
    if np.sum(['M3' in i for i in statistics])>=1:
        pos=np.where(['M3' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",1)[1])
            aux=Datos.resample(str(h) + t , how='sum');
            statististics_values_real.append(sp.stats.moment(aux, moment=3, nan_policy='omit'))

    return statististics_values_real



class evaluateInd_PSO(object):
    def __init__(self, vector, weights, IF, process, statistics):###__init__Siembre que ejecutes la clase evaluateInd_PSO va a haer ese self.v
        self.v = vector
        self.w = weights
        self.I_F = IF
        self.p = process
        self.s = statistics
    def __call__(self, ind):###__call__Si haces evaluateInd_PSO() te hace el class
        ##Defino variables
        landa=list(); mu_c=list(); eta=list(); xi=list(); betha=list(); alpha=list(); alpha_p=list()
        if self.p=='normal'and self.I_F=='E':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1); 
        elif self.p=='normal'and self.I_F=='W':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(1);  
        elif self.p=='normal'and self.I_F=='G':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(0.5); 
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(0.5); 
        elif self.p=='cells' and self.I_F=='E':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(ind[5]); 
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[6]); xi.append(ind[7]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1-ind[5]);
        elif self.p=='cells'and self.I_F=='W':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(ind[5]); 
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[6]); xi.append(ind[7]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1-ind[5]);
        elif self.p=='cells'and self.I_F=='G':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(ind[6]); 
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[7]); xi.append(ind[8]); betha.append(ind[4]); alpha.append(9); alpha_p.append(1-ind[6]);
        elif self.p=='storms'and self.I_F=='E':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1); 
            landa.append(ind[5]); mu_c.append(ind[6]); eta.append(ind[7]); xi.append(ind[8]); betha.append(ind[9]); alpha.append(1); alpha_p.append(1); 
        elif self.p=='storms'and self.I_F=='W':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(1); 
            landa.append(ind[6]); mu_c.append(ind[7]); eta.append(ind[8]); xi.append(ind[9]); betha.append(ind[10]); alpha.append(ind[11]); alpha_p.append(1); 
        elif self.p=='storms'and self.I_F=='G':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(1); 
            landa.append(ind[6]); mu_c.append(ind[7]); eta.append(ind[8]); xi.append(ind[9]); betha.append(ind[10]); alpha.append(ind[11]); alpha_p.append(1); 
            
        ##Media
        if np.sum(['mean' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            a = NSRP_mean(self.I_F,1, landa, mu_c, eta, xi, alpha, alpha_p)
            #a = self.w[0]*abs(a-self.v[0])/self.v[0]; d_e={'e0':a}
            a = self.w[0]*((1-(self.v[0]/a))**2 + (1-(a/self.v[0]))**2); d_e={'e0':a}  
            ##Eq (1) --> Regionalised spatiotemporal rainfall and temperature models for flood studies in the
            ## Basque Country, Spain
            
        ##Varianza
        if np.sum(['var' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['var' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_covariance(self.I_F,h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                #a = self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
                
        ##Autocorrelacion
        if np.sum(['autocorr' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['autocorr' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                l=int(self.s[ii].split("_",3)[1])
                h=int(self.s[ii].split("_",3)[2])
                a=NSRP_covariance(self.I_F,h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                    /NSRP_covariance(self.I_F,h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                #a=self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
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
                a = NSRP_momentos_orden_3(self.I_F,'Poisson',h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                #a=self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        
        for k in d_e.keys():
            if isnan(d_e[k]): d_e[k] = 10000

        errores=list(); errores.append([d_e['e' + str(i)] for i in range(len(d_e.keys()))])
        
        return np.sum(errores)

    def totalError(self, ind):
        ##Defino variables
        landa=list(); mu_c=list(); eta=list(); xi=list(); betha=list(); alpha=list(); alpha_p=list()
        if self.p=='normal'and self.I_F=='E':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1); 
        elif self.p=='normal'and self.I_F=='W':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(1);  
        elif self.p=='normal'and self.I_F=='G':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(0.5); 
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(0.5); 
        elif self.p=='cells' and self.I_F=='E':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(ind[5]); 
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[6]); xi.append(ind[7]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1-ind[5]);
        elif self.p=='cells'and self.I_F=='W':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(ind[5]); 
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[6]); xi.append(ind[7]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1-ind[5]);
        elif self.p=='cells'and self.I_F=='G':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(ind[6]); 
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[7]); xi.append(ind[8]); betha.append(ind[4]); alpha.append(9); alpha_p.append(1-ind[6]);
        elif self.p=='storms'and self.I_F=='E':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1); 
            landa.append(ind[5]); mu_c.append(ind[6]); eta.append(ind[7]); xi.append(ind[8]); betha.append(ind[9]); alpha.append(1); alpha_p.append(1); 
        elif self.p=='storms'and self.I_F=='W':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(1); 
            landa.append(ind[6]); mu_c.append(ind[7]); eta.append(ind[8]); xi.append(ind[9]); betha.append(ind[10]); alpha.append(ind[11]); alpha_p.append(1); 
        elif self.p=='storms'and self.I_F=='G':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(1); 
            landa.append(ind[6]); mu_c.append(ind[7]); eta.append(ind[8]); xi.append(ind[9]); betha.append(ind[10]); alpha.append(ind[11]); alpha_p.append(1);  
        
        ##Media
        if np.sum(['mean' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            a = NSRP_mean(self.I_F,1, landa, mu_c, eta, xi, alpha, alpha_p)
            #a = self.w[0]*abs(a-self.v[0])/self.v[0]; d_e={'e0':a}
            a = self.w[0]*((1-(self.v[0]/a))**2 + (1-(a/self.v[0]))**2); d_e={'e0':a}
            
        ##Varianza
        if np.sum(['var' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['var' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_covariance(self.I_F,h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                #a = self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
                
        ##Autocorrelacion
        if np.sum(['autocorr' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['autocorr' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                l=int(self.s[ii].split("_",3)[1])
                h=int(self.s[ii].split("_",3)[2])
                a=NSRP_covariance(self.I_F,h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                    /NSRP_covariance(self.I_F,h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
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
                a=NSRP_momentos_orden_3(self.I_F,'Poisson',h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                #a = self.w[ii]*abs(a-self.v[ii])/self.v[ii]; d_e['e' + str(ii)]=a
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        
        for k in d_e.keys():
            if isnan(d_e[k]): d_e[k] = 10000
                
        errores=list(); errores.append([d_e['e' + str(i)] for i in range(len(d_e.keys()))])
        
        return np.sum(errores)
        
    def compute(self, ind):
        ##Defino variables
        landa=list(); mu_c=list(); eta=list(); xi=list(); betha=list(); alpha=list(); alpha_p=list()
        if self.p=='normal'and self.I_F=='E':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1); 
        elif self.p=='normal'and self.I_F=='W':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(0.5); 
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(0.5); 
        elif self.p=='normal'and self.I_F=='G':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(0.5); 
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(0.5); 
        elif self.p=='cells' and self.I_F=='E':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(ind[5]); 
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[6]); xi.append(ind[7]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1-ind[5]);
        elif self.p=='cells'and self.I_F=='W':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(ind[5]); 
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[6]); xi.append(ind[7]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1-ind[5]);
        elif self.p=='cells'and self.I_F=='G':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(ind[6]); 
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[7]); xi.append(ind[8]); betha.append(ind[4]); alpha.append(9); alpha_p.append(1-ind[6]);
        elif self.p=='storms'and self.I_F=='E':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1); 
            landa.append(ind[5]); mu_c.append(ind[6]); eta.append(ind[7]); xi.append(ind[8]); betha.append(ind[9]); alpha.append(1); alpha_p.append(1); 
        elif self.p=='storms'and self.I_F=='W':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(1); 
            landa.append(ind[6]); mu_c.append(ind[7]); eta.append(ind[8]); xi.append(ind[9]); betha.append(ind[10]); alpha.append(ind[11]); alpha_p.append(1); 
        elif self.p=='storms'and self.I_F=='G':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(ind[5]); alpha_p.append(1); 
            landa.append(ind[6]); mu_c.append(ind[7]); eta.append(ind[8]); xi.append(ind[9]); betha.append(ind[10]); alpha.append(ind[11]); alpha_p.append(1); 
        
        ##Media
        if np.sum(['mean' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            a = NSRP_mean(self.I_F,1, landa, mu_c, eta, xi, alpha, alpha_p)
            v={'v0':a}
        ##Varianza
        if np.sum(['var' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['var' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_covariance(self.I_F,h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                v['v' + str(ii)]=a
        ##Autocorrelacion
        if np.sum(['autocorr' in i for i in self.s])>=1:#Quiere decir que hay que ajustar ese estadístico
            pos=np.where(['autocorr' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                l=int(self.s[ii].split("_",3)[1])
                h=int(self.s[ii].split("_",3)[2])
                a=NSRP_covariance(self.I_F,h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
                    /NSRP_covariance(self.I_F,h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
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
                a=NSRP_momentos_orden_3(self.I_F,'Poisson',h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                v['v' + str(ii)]=a
            
        valores=list(); valores.append([v['v' + str(i)] for i in range(len(v.keys()))])
        
        return valores


def NSRP_simulation1(Parametros_meses, year_ini, year_fin, temporal_resolution, I_F, process, Seasonality):

    
    landa=Parametros_meses[Parametros_meses.index=='landa'].values[0];Storm_origin=1/landa##Poisson process
    ipsilon=Parametros_meses[Parametros_meses.index=='ipsilon'].values[0];Number_cells_per_storm=ipsilon-1##Random poisson mean ¿ipsilon-1?
    if process=='normal':
        eta=Parametros_meses[Parametros_meses.index=='eta'].values[0]; Duracion_cell=1/eta##exponencial random
        xi=Parametros_meses[Parametros_meses.index=='xi'].values[0]; Intensidad_cell=1/xi##exponencia random###ojo
        if I_F=='W' or I_F=='G':
            alpha=Parametros_meses[Parametros_meses.index=='alpha'].values[0];
    elif process=='cells':    
        eta1=Parametros_meses[Parametros_meses.index=='eta1'].values[0]; Duracion_cell1=1/eta1##exponencial random
        xi1=Parametros_meses[Parametros_meses.index=='xi1'].values[0]; Intensidad_cell1=1/xi1##exponencia random###ojo
        eta2=Parametros_meses[Parametros_meses.index=='eta2'].values[0]; Duracion_cell2=1/eta2##exponencial random
        xi2=Parametros_meses[Parametros_meses.index=='xi2'].values[0]; Intensidad_cell2=1/xi2##exponencia random###ojo
        alpha_p1=Parametros_meses[Parametros_meses.index=='alpha_p1'].values[0]
        if I_F=='W' or I_F=='G':
            alpha1=Parametros_meses[Parametros_meses.index=='alpha1'].values[0];
            alpha2=Parametros_meses[Parametros_meses.index=='alpha2'].values[0];

    betha=Parametros_meses[Parametros_meses.index=='betha'].values[0]; Dist_cell_origin=1/betha##exponencial random


    time_d=pd.date_range(start=str((year_ini))+'-01-01',  end=str((year_fin))+'-01-01', freq='D')
    time_h=pd.date_range(start=str((year_ini))+'-01-01',  end=str((year_fin))+'-01-01', freq='h')
    time_min=pd.date_range(start=str((year_ini))+'-01-01',  end=str((year_fin))+'-01-01', freq='min')

    Dataframe_simulacion_unido_day_estaciones=pd.DataFrame(index=time_d)
    Dataframe_simulacion_unido_hour_estaciones=pd.DataFrame(index=time_h)
    Dataframe_simulacion_unido_min_estaciones=pd.DataFrame(index=time_min)
    Dataframe_simulacion_unido_day_estaciones['Rain']=np.ones(len(Dataframe_simulacion_unido_day_estaciones))*0
    Dataframe_simulacion_unido_hour_estaciones['Rain']=np.ones(len(Dataframe_simulacion_unido_hour_estaciones))*0
    tt=pd.date_range(start=datetime.datetime(year_ini,1, 1),end=datetime.datetime(year_fin,1, 1), freq='min')##Si esto lo pones en
    ##en segundos abajo hay que dividir entre 3600 y no 60
    #tt_ordinal=[datetime2matlabdn(ii) for i, ii in enumerate(tt)]##Solo hay que hacerlo una vez
    #tt_ordinal=list(np.loadtxt('tt_ordinal_10years.txt'))
    
    Intensidad_cells_total=list()
    Duracion_cells_total=list()

    for monthss in Seasonality:
        n_days=np.sum(np.in1d(Dataframe_simulacion_unido_day_estaciones.index.month, monthss))#calculo el numero de
        #dias que hay en los meses i
        position_day=np.in1d(Dataframe_simulacion_unido_day_estaciones.index.month, monthss)
        position_hour=np.in1d(Dataframe_simulacion_unido_hour_estaciones.index.month, monthss)

        Dataframe_simulacion_day_aux=pd.date_range('1/1/1800', periods=n_days, freq='D')
        Dataframe_simulacion_hour_aux=pd.date_range('1/1/1800', periods=n_days, freq='h')
        Dataframe_simulacion_min_aux=pd.date_range('1/1/1800', periods=n_days, freq='min')
        
        time_star=Dataframe_simulacion_day_aux[0]
        time_end=Dataframe_simulacion_day_aux[-1]
        time_lapso=time_star
        time_storm_origins=list()
        n=0
        
        if np.size(monthss)==1:
            monthss=[monthss]
        
        while time_lapso < time_end:
            s = np.random.exponential(Storm_origin[monthss[0]-1], 1)

            if temporal_resolution=='d':
                time_lapso=time_lapso +  datetime.timedelta(days=s[0])##CAMBIAR dependeindo si estas en h o d
            elif temporal_resolution=='h':
                time_lapso=time_lapso +  datetime.timedelta(hours=s[0])##CAMBIAR dependeindo si estas en h o d

            time_storm_origins.append(time_lapso)
            n=n+1


        n_storms=len(time_storm_origins)
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
        for i in range(n_storms):##numero de rain cells
            time1=time_storm_origins[i] #ojo horas!
            Distancia_horas_cell_sim=np.random.exponential(scale=Dist_cell_origin[monthss[0]-1], size=Number_cell_per_storm[i])

            if process=='normal':
                Duracion_horas_cell_sim=np.random.exponential(scale=Duracion_cell[monthss[0]-1], size=Number_cell_per_storm[i])
                if I_F=='E':
                    Intensidad_cell_sim=np.random.exponential(scale=Intensidad_cell[monthss[0]-1], size=Number_cell_per_storm[i])
                    #Intensidad_cell_sim=np.random.weibull(1, size=Number_cell_per_storm[i])*Intensidad_cell[monthss[0]-1]

                elif I_F=='W':
                    Intensidad_cell_sim=np.random.weibull(alpha[monthss[0]-1], size=Number_cell_per_storm[i])*Intensidad_cell[monthss[0]-1]
                elif I_F=='G':
                    Intensidad_cell_sim=np.random.gamma(alpha[monthss[0]-1], Intensidad_cell[monthss[0]-1], size=Number_cell_per_storm[i])

            for j in range(Number_cell_per_storm[i][0]):#Nuevo

                if temporal_resolution=='d':
                    time1_cell=time1 + datetime.timedelta(days=Distancia_horas_cell_sim[j])##CAMBIAR dependeindo si estas en h o d     
                    time2_cell=time1_cell + datetime.timedelta(days=Duracion_horas_cell_sim[j])##CAMBIAR dependeindo si estas en h o d

                elif temporal_resolution=='h':
                    time1_cell=time1 + datetime.timedelta(hours=Distancia_horas_cell_sim[j])##CAMBIAR dependeindo si estas en h o d
                    time2_cell=time1_cell + datetime.timedelta(hours=Duracion_horas_cell_sim[j])##CAMBIAR dependeindo si estas en h o d

                Intensidad=Intensidad_cell_sim[j]
                time_ini_cells.append(time1_cell)
                time_fin_cells.append(time2_cell)
                Intensidad_cells.append(Intensidad)
                Duracion_horas_cells.append(Duracion_horas_cell_sim[j])
        ################################################################################
        time_ini_cells=np.array(time_ini_cells)
        time_fin_cells=np.array(time_fin_cells)
        Intensidad_cells=np.array(Intensidad_cells)
        #x_cells=np.array(x_cells); y_cells=np.array(y_cells)
        #radio_cells=np.array(radio_cells)
        Duracion_horas_cells=np.array(Duracion_horas_cells)
        #################################################################################
        ############################################################################
        tt=pd.date_range(start=Dataframe_simulacion_day_aux[0],end=Dataframe_simulacion_day_aux[-1], freq='min')##Si esto lo pones en
        ##en segundos abajo hay que dividir entre 3600 y no 60
        tt_ordinal=datetime2matlabdnJavi(tt, 'minutes')
        #tt_ordinal=[datetime2matlabdn(ii) for i, ii in enumerate(tt)]##Solo hay que hacerlo una vez
        ############################################################################
        Anhos=list()
        for i, ii in enumerate(time_fin_cells):
            aux=ii; year_aux=aux.year
            Anhos.append(year_aux)
        ##Quito las cenldas que dentro de una tormenta caen fuera del año límite
        Dentro_fechass=np.where(np.array(Anhos)<=Dataframe_simulacion_day_aux.year[-1])
        time_ini_cellss=time_ini_cells[Dentro_fechass[0]]
        time_fin_cellss=time_fin_cells[Dentro_fechass[0]]
        Intensidad_cellss=Intensidad_cells[Dentro_fechass[0]]; Intensidad_cells_total.append(Intensidad_cellss)
        Duracion_horas_cellss=Duracion_horas_cells[Dentro_fechass[0]]; Duracion_cells_total.append(Duracion_horas_cellss)
        ##Acumulo cada celda de lluvia en minutos para luego agruparlo en horas
        zeros=np.zeros((len(Dataframe_simulacion_min_aux), 1))
        aux_t=np.zeros((len(Dataframe_simulacion_min_aux), 1))
        t_ini = np.hstack([time_ini_cellss, time_fin_cellss])
        if process=='cells':
            i_ini = np.hstack([Intensidad_cellss.T, -Intensidad_cellss.T])
            i_ini = i_ini[0]
        else: 
            i_ini = np.hstack([Intensidad_cellss, -Intensidad_cellss])
        orden = np.argsort(t_ini)
        t = t_ini[orden]
        i = np.cumsum(i_ini[orden])
        i[i<0] = 0
        rain=i.copy()
        t_ordinal=[datetime2matlabdn(ii) for i, ii in enumerate(t)]

        rainfall = interp1d(t_ordinal, rain, kind="zero", bounds_error=False, fill_value=0.)
        rr = rainfall(tt_ordinal)
        Date=pd.DataFrame(index=tt)
        if temporal_resolution=='d':
            Date['Rain']=rr/(60*24)
        elif temporal_resolution=='h':
            Date['Rain']=rr/(60)###calculo todo en horas y lo estoy interpolando a minutos por eso lo divido entre 60. Si interpolase
        ##a segundos habria que utilizar poner otro 60. Si las tormentas son muy pequeñas pierdo algo de información haciendolo
        ##en minutos...El problma es que en segundos tarda una eternidad
        #Date['Rain']=rr/3600###calculo todo en horas y lo estoy interpolando a minutos por eso lo divido entre 60. Si interpolase
        ##a segundos habria que utilizar poner otro 60
        Dataframe_simulacion_unido_day_estaciones['Rain'].\
        iloc[np.where(position_day)[0][0:len(Date.resample('D', how='sum').values)]]=\
        Date.resample('D', how='sum').values
        Dataframe_simulacion_unido_hour_estaciones['Rain'].\
        iloc[np.where(position_hour)[0][0:len(Date.resample('h', how='sum').values)]]=\
        Date.resample('h', how='sum').values
        del Date, rr, rainfall, t_ordinal, rain, i, t, orden, t_ini, i_ini, aux_t, zeros
    return Dataframe_simulacion_unido_hour_estaciones, Dataframe_simulacion_unido_day_estaciones,\
           np.hstack(Intensidad_cells_total), np.hstack(Duracion_cells_total)



