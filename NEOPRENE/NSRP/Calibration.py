'''
Library containing classes for calibrating model parameters.

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''

from NEOPRENE.NSRP.MathematicalPropertiesNSRP import *
from NEOPRENE.NSRP.utils import *
from NEOPRENE.NSRP.libs_NSRP import *
import time
import yaml
import ast
from NEOPRENE.NSRP.parameters_calibration import parameters_calibration

class Calibration(object):
    """This function allows you to configure the point model to calibrate from user-defined hyperparameters.
    
    Parameters
    ----------
    hiperparams : Object. 
        Object with calibration hyperparameters.
      
    """
    def __init__(self,hiperparams):
        self.hiperparams = hiperparams
        
    def __call__(self, statistics, verbose=False):
        
        """Starting from the general hyperparameters of the model and different real statistics, the calibration of the point model begins.
        
        Parameters
        ----------
        statistics: Object. 
            Object with statistics.
          
        """
        
        
    
        self.statististics_dataframe   = statistics.statististics_dataframe
   
        ##Limits are obtained at the indicated time resolution
        if self.hiperparams.temporal_resolution=='d':
            t=24
        else:
            t=1
            
        lim = np.array([
            [(1/self.hiperparams.time_between_storms[1])*t, (1/self.hiperparams.time_between_storms[0])*t],
            [self.hiperparams.number_storm_cells[0], self.hiperparams.number_storm_cells[1]],
            [(1/self.hiperparams.cell_duration[1])*t, (1/self.hiperparams.cell_duration[0])*t],
            [(1/self.hiperparams.cell_intensity[1])*t, (1/self.hiperparams.cell_intensity[0])*t], 
            [(1/self.hiperparams.storm_cell_displacement[1])*t, (1/self.hiperparams.storm_cell_displacement[0])*t]])
            
        statististics_fit_dataframe=pd.DataFrame(index=self.hiperparams.statistics_name,columns=self.hiperparams.Seasonality_str)

        if self.hiperparams.process=='normal':
            param_s=['landa', 'ipsilon', 'eta', 'xi', 'betha'] 
            
            # We create a dataframe to store the adjusted intensity function parameters.
            Dataframe_params=pd.DataFrame(index=param_s,columns=self.hiperparams.Seasonality_str); n_p=5;
            limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                              [lim[4,0],lim[4,1]]])

        elif self.hiperparams.process=='storms':
            param_s=['landa1', 'ipsilon1', 'eta1', 'xi1', 'betha1', 'landa2', 'ipsilon2', 'eta2', 'xi2', 'betha2']
            
            # We create a dataframe to store the adjusted intensity function parameters.
            Dataframe_params=pd.DataFrame(index=param_s,columns=self.hiperparams.Seasonality_str); n_p=10
            limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                              [lim[4,0],lim[4,1]], [lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],\
                              [lim[3,0],lim[3,1]],[lim[4,0],lim[4,1]]]) 
                              
        ## We start the adjustment for each with seasonality.
        print('')
        print('')
        print('#'*80)
        print('Adjustment of parameters using the Particle Swarm Optimization (PSO)')
        print('')
        print('')
        Error_dataframe=pd.DataFrame(index=self.hiperparams.Seasonality_str)
        Error_dataframe['Total Error']=np.nan
        
        for pri, prii in enumerate(self.hiperparams.Seasonality):
            evaluator = evaluateInd_PSO(self.statististics_dataframe.loc[:,str(prii)].values, self.hiperparams.weights, self.hiperparams.process, self.hiperparams.statistics_name)
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
                statististics_fit_dataframe.loc[:,str(prii)]= sta_fit[0]
                Fitted_parameters=mySwarm.getBest()[1]
                Error_dataframe.loc[str(prii)]=evaluator.totalError(bp)

                total_error=evaluator.totalError(bp)

            if n_i>1:
                # If in the first initialisation it does not find the best parameters, it does it again and keeps the simulation with the lowest error.
                Error_dataframe.loc[str(prii)]=np.min(error_initializations)
                sta_fit=evaluator.compute(parameters_initializations[np.argmin(error_initializations)]);
                statististics_fit_dataframe.loc[:,str(prii)]= sta_fit[0]
                Fitted_parameters=parameters_initializations[np.argmin(error_initializations)]


            if self.hiperparams.process=='normal':

                param_v=list(); 
                param_v.append(Fitted_parameters[0]);
                param_v.append(Fitted_parameters[1]); 
                param_v.append(Fitted_parameters[2]); 
                param_v.append(Fitted_parameters[3]); 
                param_v.append(Fitted_parameters[4]); 
                Dataframe_params.loc[:,str(prii)]=param_v



            elif self.hiperparams.process=='storms':

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
                Dataframe_params.loc[:,str(prii)]=param_v

        resuls = parameters_calibration(Dataframe_params,statististics_fit_dataframe,self.statististics_dataframe, Error_dataframe)
        
        
        return resuls
    
