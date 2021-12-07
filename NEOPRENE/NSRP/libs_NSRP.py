'''
Library containing the necessary functions to simulate the Neyman-Scott process. 
Rectangular Pulses (NSRP).

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''

'''Functions of the NSRP mode'''

from NEOPRENE.NSRP.MathematicalPropertiesNSRP import *
from NEOPRENE.NSRP.utils import *
import NEOPRENE.NSRP.PSOMdJ as pso
import numpy as np
import pandas as pd
import scipy as sp
import math as mt
import datetime
from scipy.interpolate import interp1d
from datetime import date
from datetime import timedelta
from scipy import stats


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

    if np.sum(['mean' in i for i in statistics])>=1:
        statististics_values_real.append(np.nanmean(Data))
    if np.sum(['var' in i for i in statistics])>=1:
        pos=np.where(['var' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",1)[1])
            aux=Data.resample(str(h) + t).agg(pd.Series.sum, min_count=1); 
            statististics_values_real.append(np.nanvar(aux))
    if np.sum(['autocorr' in i for i in statistics])>=1:
        pos=np.where(['autocorr' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            l=int(statistics[ii].split("_",3)[1])
            h=int(statistics[ii].split("_",3)[2])
            aux=Data.resample(str(h) + t).agg(pd.Series.sum, min_count=1); 
            Autocorrelation_aux=aux.autocorr(lag=l) 
            if np.size(Autocorrelation_aux)>1: Autocorrelation_aux=Autocorrelation_aux[0] 
            statististics_values_real.append(Autocorrelation_aux)
    if np.sum(['fih' in i for i in statistics])>=1:
        pos=np.where(['fih' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",1)[1])
            statististics_values_real.append(fi_h(Data, h))
    if np.sum(['fiWW' in i for i in statistics])>=1:
        pos=np.where(['fiWW' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",1)[1])
            statististics_values_real.append(fi_WW(Data, h))
    if np.sum(['fiDD' in i for i in statistics])>=1:
        pos=np.where(['fiDD' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",1)[1])
            statististics_values_real.append(fi_DD(Data, h))
    if np.sum(['M3' in i for i in statistics])>=1:
        pos=np.where(['M3' in i for i in statistics]); pos=pos[0]
        for i, ii in enumerate(pos):
            h=int(statistics[ii].split("_",1)[1])
            aux=Data.resample(str(h) + t ).agg(pd.Series.sum, min_count=1);
            statististics_values_real.append(sp.stats.moment(aux, moment=3, nan_policy='omit'))
    
    return statististics_values_real



class evaluateInd_PSO(object):
    def __init__(self, vector, weights, process, statistics):
        self.v = vector
        self.w = weights
        self.p = process
        self.s = statistics
    def __call__(self, ind):
    
        landa=list(); mu_c=list(); eta=list(); xi=list(); betha=list(); alpha=list(); alpha_p=list()
        if self.p=='normal':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1); 

        elif self.p=='storms':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1); 
            landa.append(ind[5]); mu_c.append(ind[6]); eta.append(ind[7]); xi.append(ind[8]); betha.append(ind[9]); alpha.append(1); alpha_p.append(1); 
            
        ##Mean
        if np.sum(['mean' in i for i in self.s])>=1:
            a = NSRP_mean(1, landa, mu_c, eta, xi, alpha, alpha_p)
            a = self.w[0]*((1-(self.v[0]/a))**2 + (1-(a/self.v[0]))**2); d_e={'e0':a}  

            
        ##Variance
        if np.sum(['var' in i for i in self.s])>=1:
            pos=np.where(['var' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        
        ##Autocorrlation
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
                a = NSRP_moments_order_3('Poisson',h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)

                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        
        for k in d_e.keys():
            if isnan(d_e[k]): d_e[k] = 10000

        error=list(); error.append([d_e['e' + str(i)] for i in range(len(d_e.keys()))])
        
        return np.sum(error)

    def totalError(self, ind):

        landa=list(); mu_c=list(); eta=list(); xi=list(); betha=list(); alpha=list(); alpha_p=list()
        if self.p=='normal':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1); 
        elif self.p=='storms':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1); 
            landa.append(ind[5]); mu_c.append(ind[6]); eta.append(ind[7]); xi.append(ind[8]); betha.append(ind[9]); alpha.append(1); alpha_p.append(1); 
        
        ##Mean
        if np.sum(['mean' in i for i in self.s])>=1:
            a = NSRP_mean(1, landa, mu_c, eta, xi, alpha, alpha_p)
            a = self.w[0]*((1-(self.v[0]/a))**2 + (1-(a/self.v[0]))**2); d_e={'e0':a}
            
        ##Variance
        if np.sum(['var' in i for i in self.s])>=1:
            pos=np.where(['var' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
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
                a=NSRP_moments_order_3('Poisson',h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                a = self.w[ii]*((1-(self.v[ii]/a))**2 + (1-(a/self.v[ii]))**2); d_e['e' + str(ii)]=a
        
        for k in d_e.keys():
            if isnan(d_e[k]): d_e[k] = 10000
                
        error=list(); error.append([d_e['e' + str(i)] for i in range(len(d_e.keys()))])
        
        return np.sum(error)
        
    def compute(self, ind):
        landa=list(); mu_c=list(); eta=list(); xi=list(); betha=list(); alpha=list(); alpha_p=list()
        if self.p=='normal':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1); 
        elif self.p=='storms':
            landa.append(ind[0]); mu_c.append(ind[1]); eta.append(ind[2]); xi.append(ind[3]); betha.append(ind[4]); alpha.append(1); alpha_p.append(1); 
            landa.append(ind[5]); mu_c.append(ind[6]); eta.append(ind[7]); xi.append(ind[8]); betha.append(ind[9]); alpha.append(1); alpha_p.append(1); 
        
        ##Mean
        if np.sum(['mean' in i for i in self.s])>=1:
            a = NSRP_mean(1, landa, mu_c, eta, xi, alpha, alpha_p)
            v={'v0':a}
        ##Variance
        if np.sum(['var' in i for i in self.s])>=1:
            pos=np.where(['var' in i for i in self.s]); pos=pos[0]
            for i, ii in enumerate(pos):
                h=int(self.s[ii].split("_",1)[1])
                a = NSRP_covariance(h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
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
                a=NSRP_moments_order_3('Poisson',h,0, landa, mu_c, eta, xi, betha, alpha, alpha_p)
                v['v' + str(ii)]=a
            
        values=list(); values.append([v['v' + str(i)] for i in range(len(v.keys()))])
        
        return values
	
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
            i_ini = np.hstack([Intensity_cellss.T, -Intensity_cellss.T])
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






