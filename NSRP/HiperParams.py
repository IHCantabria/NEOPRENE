'''
Library containing classes for reading yaml configuration files for calibration and simulation.

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''



import yaml
import ast

class Calibration(object):
    """This function allows to obtain the hyperparameters to start the calibration of the point model and to check if all parameters are correct.
    
    Parameters
    ----------
    file_yml : String. 
        Directory where the Yml file with the configuration to perform the model calibration is located.
      
    """
    def __init__(self,file_yml):
        with open(file_yml) as file:
            hiper_params = yaml.load(file, Loader=yaml.FullLoader)
        self.Seasonality_type        = hiper_params['Seasonality_type']
        self.Seasonality_user        = hiper_params['Seasonality_user']
        self.temporal_resolution     = hiper_params['temporal_resolution']
        self.process                 = hiper_params['process']
        self.statistics_name         = hiper_params['statistics_name']
        self.weights                 = hiper_params['weights']
        self.number_iterations       = hiper_params['number_iterations']
        self.number_bees             = hiper_params['number_bees']
        self.number_initializations  = hiper_params['number_initializations']
        self.time_between_storms     = hiper_params['time_between_storms']
        self.number_storm_cells      = hiper_params['number_storm_cells']
        self.cell_duration           = hiper_params['cell_duration']
        self.cell_intensity          = hiper_params['cell_intensity']
        self.storm_cell_displacement = hiper_params['storm_cell_displacement']
        
        
        if self.Seasonality_type == 'annual':
            Seasonality=list()
            Seasonality.append((1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ,11, 12))
            
            Seasonality_str = list()
            for sr in Seasonality:
                Seasonality_str.append(str(sr))
            
        elif self.Seasonality_type == 'seasonal':
            Seasonality=list()
            Seasonality.append((1, 2, 3))
            Seasonality.append((4, 5, 6))
            Seasonality.append((7, 8, 9))
            Seasonality.append((10, 11, 12))
            
            Seasonality_str = list()
            for sr in Seasonality:
                Seasonality_str.append(str(sr))
            
        elif self.Seasonality_type == 'monthly':
            Seasonality=list()
            Seasonality.append((1))
            Seasonality.append((2))
            Seasonality.append((3))
            Seasonality.append((4))
            Seasonality.append((5))
            Seasonality.append((6))
            Seasonality.append((7))
            Seasonality.append((8))
            Seasonality.append((9))
            Seasonality.append((10))
            Seasonality.append((11))
            Seasonality.append((12))
            
            Seasonality_str = list()
            for sr in Seasonality:
                Seasonality_str.append(str(sr))
                
        elif self.Seasonality_type == 'user_defined':
        
            Seasonality_str = hiper_params['Seasonality_user']
            
            Seasonality = list()
            for i in self.Seasonality_user:
                Seasonality.append(ast.literal_eval(i))
                
            if type(hiper_params['Seasonality_user'][0])!=str:
                raise Exception ('The definition of the months that make up the user-defined seasonality is incorrect, check that it is a string.')
                
        else:
            raise Exception ('The name of the seasonality type is not correct, check that the name is among those listed below: annual, monthly, seasonal, user_defined')       
        self.Seasonality_str = Seasonality_str        
        self.Seasonality = Seasonality
        
        
            
        if self.temporal_resolution != 'h' and self.temporal_resolution != 'd':
            raise Exception ("The time resolution of the point model is poorly defined. It should be 'h' or 'd'.")
        
        if self.process != 'storms' and self.process != 'normal':
            raise Exception ("The process type is wrongly defined, it should be 'storms' or 'normal'.")
            
            
        statistic_list = ['mean_h', 'var_h', 'autocorr_l_h','fih_h', 'fiWW_h', 'fiDD_h', 'M3_h']
        
        for i in self.statistics_name:
            if i.split('_')[0]== 'autocorr':
                st =  i.split('_')[0]+'_'+'l'+'_'+'h'
                if st not in statistic_list:
                    raise Exception ('The '+i+' statistic has not been correctly defined.')
            else:
                st = i.split('_')[0]+'_'+'h'
                if st not in statistic_list:
                    raise Exception ('The '+i+' statistic has not been correctly defined.')
        
            
            
            
        
        
        
class Simulation(object): 
    """This function allows to obtain the hyperparameters to start the simulation of the point model and to check if all parameters are correct.
        
        Parameters
        ----------
        file_yml : String. 
            Directory where the Yml file with the configuration to perform the model simulation is located.
    """
    def __init__(self,file_yml_sim): 
        with open(file_yml_sim) as file:
            hiper_params = yaml.load(file, Loader=yaml.FullLoader)
        self.Seasonality_type        = hiper_params['Seasonality_type']
        self.Seasonality_user        = hiper_params['Seasonality_user']
        self.statistics_name         = hiper_params['statistics_name']
        self.process                 = hiper_params['process']
        self.temporal_resolution     = hiper_params['temporal_resolution']
        self.year_ini                = hiper_params['year_ini']
        self.year_fin                = hiper_params['year_fin']
        
        if self.Seasonality_type=='annual':
            Seasonality=list()
            Seasonality.append((1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ,11, 12))
            
            Seasonality_str = list()
            for sr in Seasonality:
                Seasonality_str.append(str(sr))
            
        elif self.Seasonality_type=='seasonal':
            Seasonality=list()
            Seasonality.append((1, 2, 3))
            Seasonality.append((4, 5, 6))
            Seasonality.append((7, 8, 9))
            Seasonality.append((10, 11, 12))
            
            Seasonality_str = list()
            for sr in Seasonality:
                Seasonality_str.append(str(sr))
            
        elif self.Seasonality_type =='monthly':
            Seasonality=list()
            Seasonality.append((1))
            Seasonality.append((2))
            Seasonality.append((3))
            Seasonality.append((4))
            Seasonality.append((5))
            Seasonality.append((6))
            Seasonality.append((7))
            Seasonality.append((8))
            Seasonality.append((9))
            Seasonality.append((10))
            Seasonality.append((11))
            Seasonality.append((12))
            
            Seasonality_str = list()
            for sr in Seasonality:
                Seasonality_str.append(str(sr))
                
        elif self.Seasonality_type =='user_defined':
        
            Seasonality_str = hiper_params['Seasonality_user']
            
            Seasonality = list()
            for i in self.Seasonality_user:
                Seasonality.append(ast.literal_eval(i))
                
            if type(hiper_params['Seasonality_user'][0])!=str:
                raise Exception ('The definition of the months that make up the user-defined seasonality is incorrect, check that it is a string.')
                
        else:
            raise Exception ('The name of the seasonality type is not correct, check that the name is among those listed below: annual, monthly, seasonal, user_defined')
                
        self.Seasonality_str = Seasonality_str        
        self.Seasonality = Seasonality