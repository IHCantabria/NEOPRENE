'''
Library containing classes for calibrating model parameters.

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''

from NEOPRENE.STNSRP.MathematicalPropertiesSTNSRP import *
from NEOPRENE.STNSRP.utils import *
from NEOPRENE.STNSRP.libs_STNSRP import *
import time
import yaml
import ast
from NEOPRENE.STNSRP.parameters_calibration import parameters_calibration
import sys

class Calibration(object):
    """This function allows you to configure the STNSRPM to calibrate from user-defined hyperparameters.
    
    Parameters
    ----------
    hiperparams : Object. 
        Object with calibration hyperparameters.
      
    """
    def __init__(self,hiperparams):
        self.hiperparams  = hiperparams
    def __call__(self, statistics, time_series, verbose=False):
        """Starting from the general hyperparameters of the model and different real statistics, the calibration of the point model begins.
        
        Parameters
        ----------
        statistics: Object. 
            Object with statistics.
          
        """
        #self.Input_Series = Input_Series.copy()
        
    
        self.statistics_dataframe     = statistics.statistics_dataframe.copy()
        self.crosscorr_dataframe      = statistics.crosscorr_dataframe.copy()
        self.crosscorr_dataframe_dist = statistics.crosscorr_dist_dataframe.copy()
   
        ##Limits are obtained at the indicated time resolution
        if self.hiperparams.temporal_resolution=='d':
            t=24
        else:
            t=1
            
        if self.hiperparams.storm_radius==True:
            
            lim = np.array([
                [(1/self.hiperparams.time_between_storms[1])*t, (1/self.hiperparams.time_between_storms[0])*t],
                [self.hiperparams.number_storm_cells[0], self.hiperparams.number_storm_cells[1]],
                [(1/self.hiperparams.cell_duration[1])*t, (1/self.hiperparams.cell_duration[0])*t],
                #[(1/self.hiperparams.cell_intensity[1])*t, (1/self.hiperparams.cell_intensity[0])*t], 
                [(1/self.hiperparams.storm_cell_displacement[1])*t, (1/self.hiperparams.storm_cell_displacement[0])*t],
                [(1/self.hiperparams.cell_radius[1]), (1/self.hiperparams.cell_radius[0])],
                [(1/self.hiperparams.storm_radius_p[1]), (1/self.hiperparams.storm_radius_p[0])]])
            
        elif self.hiperparams.storm_radius==False:
            
              lim = np.array([
                [(1/self.hiperparams.time_between_storms[1])*t, (1/self.hiperparams.time_between_storms[0])*t],
                [self.hiperparams.number_storm_cells[0], self.hiperparams.number_storm_cells[1]],
                [(1/self.hiperparams.cell_duration[1])*t, (1/self.hiperparams.cell_duration[0])*t],
                #[(1/self.hiperparams.cell_intensity[1])*t, (1/self.hiperparams.cell_intensity[0])*t], 
                [(1/self.hiperparams.storm_cell_displacement[1])*t, (1/self.hiperparams.storm_cell_displacement[0])*t],
                [(1/self.hiperparams.cell_radius[1]), (1/self.hiperparams.cell_radius[0])]])
            
            
        #statististics_fit_dataframe=pd.DataFrame(index=self.hiperparams.statistics_name,columns=self.hiperparams.Seasonality_str)
        sta_len = len(self.hiperparams.statistics_name)
        nn_cross = np.sum(['cross' in i for i in self.hiperparams.statistics_name])
        cross_corr_division = 10#
        statistics_dataframe_fit=pd.DataFrame(index = self.hiperparams.statistics_name[0:sta_len-nn_cross])
        #crosscorr_dataframe_fit=self.crosscorr_dataframe
        # create empty dictionary where overwritte the fitted results
        crosscorr_dataframe_fit = {}
        for key in self.crosscorr_dataframe.keys():
            crosscorr_dataframe_season = {}
            for key_season in self.crosscorr_dataframe[key].keys():
                dataframe_aux = pd.DataFrame(index = self.crosscorr_dataframe[key][key_season].index, columns = self.crosscorr_dataframe[key][key_season].columns)
                dataframe_aux['dist'] = self.crosscorr_dataframe[key][key_season]['dist']
                crosscorr_dataframe_season[key_season] = dataframe_aux
            crosscorr_dataframe_fit[key] = crosscorr_dataframe_season

        if self.hiperparams.process=='normal':

            if self.hiperparams.storm_radius==False:

                param_s=['landa', 'ipsilon', 'eta', 'betha', 'fi_may'] 
                
                # We create a dataframe to store the adjusted intensity function parameters.
                Dataframe_params=pd.DataFrame(index=param_s); n_p=5;
                limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                                  [lim[4,0],lim[4,1]]])

            elif self.hiperparams.storm_radius==True:

                param_s=['landa', 'ipsilon', 'eta', 'betha', 'fi_may', 'fi_may_s'] 
                
                # We create a dataframe to store the adjusted intensity function parameters.
                Dataframe_params=pd.DataFrame(index=param_s); n_p=6;
                limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                                  [lim[4,0],lim[4,1]], [lim[5,0],lim[5,1]]])

        elif self.hiperparams.process=='storms':

            if self.hiperparams.storm_radius==False:

                param_s=['landa1', 'ipsilon1', 'eta1', 'betha1', 'fi_may1', 'landa2', 'ipsilon2', 'eta2', 'betha2', 'fi_may2']
                
                # We create a dataframe to store the adjusted intensity function parameters.
                Dataframe_params=pd.DataFrame(index=param_s); n_p=10
                limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                                  [lim[4,0],lim[4,1]],\
                                  [lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                                  [lim[4,0],lim[4,1]]])

            elif self.hiperparams.storm_radius==True:

                param_s=['landa1', 'ipsilon1', 'eta1', 'betha1', 'fi_may1', 'fi_may_s1', 'landa2', 'ipsilon2', 'eta2', 'betha2', 'fi_may1', 'fi_may_s2']
                
                # We create a dataframe to store the adjusted intensity function parameters.
                Dataframe_params=pd.DataFrame(index=param_s); n_p=12
                limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                                  [lim[4,0],lim[4,1]],[lim[5,0],lim[5,1]],\
                                  [lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                                  [lim[4,0],lim[4,1]],[lim[5,0],lim[5,1]]]) 
                              
        ## We start the adjustment for each with seasonality.
        print('')
        print('')
        print('#'*80)
        print('Adjustment of parameters using the Particle Swarm Optimization (PSO)')
        print('')
        print('')
        
        Error_dataframe=pd.DataFrame(index=self.hiperparams.Seasonality_str)
        Error_dataframe['Total Error']=np.nan

        statistics_observed_dataframe=pd.DataFrame(index=self.hiperparams.statistics_name)

        
        for pri, prii in enumerate(self.hiperparams.Seasonality):
            cross_corr_stationality_month={}
            for k in self.crosscorr_dataframe.keys():
                cross_corr_stationality_month[k]=self.crosscorr_dataframe[k][prii]
            statistics_observed_dataframe[prii]=self.statistics_dataframe[prii].mean(axis=1)
            #evaluator = evaluateInd_PSO(self.statististics_dataframe.loc[:,str(prii)].values, self.hiperparams.weights, self.hiperparams.process, self.hiperparams.statistics_name)
            evaluator = evaluateInd_PSO(self.statistics_dataframe[prii].mean(axis=1), self.hiperparams.weights, self.hiperparams.process, self.hiperparams.statistics_name, 
                                        self.hiperparams.storm_radius, cross_corr_stationality_month)
            if verbose:
                print('')
                print('Fitting the months = ' + str(prii))
                print('')
            error=0.5
            error_initializations=list()
            parameters_initializations=list()
            n_i=0
            total_error=10000
            while total_error>error and self.hiperparams.number_initializations>n_i:
                if verbose:
                    print('Initialization number  = ' + str(n_i))
                n_i+=1
                Steps = 0
                mySwarm = pso.Swarm(n_p, limites, self.hiperparams.number_bees, evaluator, 'min', verbose)
                
                while True:
                    if verbose:
                        print('Iteration number  = ' + str(Steps))
                    mySwarm.step()
                    Steps += 1
                    (b, bp) = mySwarm.getBest()
                    if verbose:
                        print('Total error = ' + str(evaluator.totalError(bp)))
                        if verbose:
                            if total_error <= error:
                                print('Stopped calibration since total error is equal or lower than error bound (' + str(error) + ')') 

                    if (evaluator.totalError(bp) < error): break
                    elif Steps > self.hiperparams.number_iterations: ###1000
                        break
                error_initializations.append(evaluator.totalError(bp))	
                parameters_initializations.append(bp)
                sta_fit=evaluator.compute(bp);
                #statististics_fit_dataframe.loc[:,str(prii)]= sta_fit[0]
                statistics_dataframe_fit.loc[:,str(prii)]=sta_fit[0][0:sta_len-nn_cross]
                for nn_ccc, st_cross in enumerate(self.hiperparams.statistics_name[sta_len-nn_cross:]):
                    crosscorr_dataframe_fit[st_cross][prii]['cross_corr'] = sta_fit[0][(sta_len-nn_cross)+(nn_ccc*cross_corr_division):(sta_len-nn_cross)+((nn_ccc+1)*cross_corr_division)]
                Fitted_parameters=mySwarm.getBest()[1]
                Error_dataframe.loc[str(prii)]=evaluator.totalError(bp)

                total_error=evaluator.totalError(bp)

            if n_i>1:
                # If in the first initialisation it does not find the best parameters, it does it again and keeps the simulation with the lowest error.
                Error_dataframe.loc[str(prii)]=np.min(error_initializations)
                sta_fit=evaluator.compute(parameters_initializations[np.argmin(error_initializations)]);
                #statististics_fit_dataframe.loc[:,str(prii)]= sta_fit[0]
                statistics_dataframe_fit.loc[:,str(prii)]=sta_fit[0][0:sta_len-nn_cross]
                for nn_ccc, st_cross in enumerate(self.hiperparams.statistics_name[sta_len:]):
                    crosscorr_dataframe_fit[st_cross][prii]['cross_corr'] = sta_fit[0][(sta_len-nn_cross)+(nn_ccc*cross_corr_division):(sta_len-nn_cross)+((nn_ccc+1)*cross_corr_division)]
                Fitted_parameters=parameters_initializations[np.argmin(error_initializations)]



            if self.hiperparams.process=='normal' and self.hiperparams.storm_radius==False:

                param_v=list(); 
                param_v.append(Fitted_parameters[0]);
                param_v.append(Fitted_parameters[1]); 
                param_v.append(Fitted_parameters[2]); 
                param_v.append(Fitted_parameters[3]); 
                param_v.append(Fitted_parameters[4]); 
                Dataframe_params[prii]=param_v

            if self.hiperparams.process=='normal' and self.hiperparams.storm_radius==True:

                param_v=list(); 
                param_v.append(Fitted_parameters[0]);
                param_v.append(Fitted_parameters[1]); 
                param_v.append(Fitted_parameters[2]); 
                param_v.append(Fitted_parameters[3]); 
                param_v.append(Fitted_parameters[4]); 
                param_v.append(Fitted_parameters[5])
                Dataframe_params[prii]=param_v


            elif self.hiperparams.process=='storms' and self.hiperparams.storm_radius==False:

                param_v=list();
                param_v.append(Fitted_parameters[0]);
                param_v.append(Fitted_parameters[1]); 
                param_v.append(Fitted_parameters[2]); 
                param_v.append(Fitted_parameters[3]); 
                param_v.append(Fitted_parameters[4]); 
                param_v.append(Fitted_parameters[5]); 
                param_v.append(Fitted_parameters[6]); 
                param_v.append(Fitted_parameters[7]); 
                param_v.append(Fitted_parameters[8]); 
                param_v.append(Fitted_parameters[9]);
                Dataframe_params[prii]=param_v
                
            elif self.hiperparams.process=='storms' and self.hiperparams.storm_radius==True:

                param_v=list();
                param_v.append(Fitted_parameters[0]);
                param_v.append(Fitted_parameters[1]); 
                param_v.append(Fitted_parameters[2]); 
                param_v.append(Fitted_parameters[3]); 
                param_v.append(Fitted_parameters[4]); 
                param_v.append(Fitted_parameters[5]); 
                param_v.append(Fitted_parameters[6]); 
                param_v.append(Fitted_parameters[7]); 
                param_v.append(Fitted_parameters[8]); 
                param_v.append(Fitted_parameters[9]); 
                param_v.append(Fitted_parameters[10]); 
                param_v.append(Fitted_parameters[11]); 
                Dataframe_params[prii]=param_v
        
        Dataframe_xi_months = XI_MONTHS(time_series, Dataframe_params, self.hiperparams.process) #calculate scale parameter for every gauge

        resuls = parameters_calibration(Dataframe_params,statistics_dataframe_fit, crosscorr_dataframe_fit, self.statistics_dataframe, self.crosscorr_dataframe,self.crosscorr_dataframe_dist, Error_dataframe, Dataframe_xi_months)
        
        
        
        
        return resuls
    
