'''
Library containing classes for reading yaml configuration files for calibration and simulation.

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''



import yaml
import ast
import sys
import numpy as np

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

        statistics_name_string = ['mean', 'var', 'autocorr','fih', 'fiWW', 'fiDD', 'M3']
        temporal_resolution_string = ['d', 'h']
        process_string = ['normal', 'storms']

        if self.process not in process_string:
            raise Exception ('The name of the hyperparameter <process> is not correct, check that the name is among those listed below: normal, process')
            sys.exit(1)  

        if self.temporal_resolution not in temporal_resolution_string:
            raise Exception ('The name of the hyperparameter <temporal_resolution> is not correct, check that the name is among those listed below: h, d')
            sys.exit(1) 

        if not len(self.time_between_storms) == 2:
            raise Exception ('The lenght of the hyperparameter <time_between_storms> must to be equal to 2')
            sys.exit(1) 

        if not len(self.number_storm_cells) == 2:
            raise Exception ('The lenght of the hyperparameter <number_storm_cells> must to be equal to 2')
            sys.exit(1) 

        if not len(self.cell_duration) == 2:
            raise Exception ('The lenght of the hyperparameter <cell_duration> must to be equal to 2')
            sys.exit(1) 

        if not len(self.cell_intensity) == 2:
            raise Exception ('The lenght of the hyperparameter <cell_intensity> must to be equal to 2')
            sys.exit(1) 

        if not len(self.storm_cell_displacement) == 2:
            raise Exception ('The lenght of the hyperparameter <storm_cell_displacement> must to be equal to 2')
            sys.exit(1) 

        for stn in self.statistics_name:
            if  not 'mean' in stn:
                if not '_' in stn:
                    raise Exception ('Lag (l) or aggregation level (h) is not included for the statistic <' + stn + '>. (e.g. ' + stn +'_1)')
                    sys.exit(1) 
                else:
                    if stn.split('_')[0] not in statistics_name_string:
                        raise Exception ('The name of the statistic <' + stn.split('_')[0] +  '> is not correct, check that the name is among those listed below: mean_h', 'var_h', 'autocorr_l_h','fih_h', 'fiWW_h', 'fiDD_h', 'M3_h')
                        sys.exit(1) 
                    else:
                        if stn.split('_')[0] in ['var','fih', 'fiWW', 'fiDD', 'M3']:
                            if len(stn.split('_'))==2:
                                if not str(stn.split('_')[1]).isnumeric():
                                    raise Exception ('Aggregation level (h) fot the statistic <' + stn.split('_')[0] +  '> has to be a integer (e.g. ' + stn +'_1)')
                                    sys.exit(1)
                            else:
                                raise Exception ('The statistic <' + stn.split('_')[0] +  '> only accept  the aggregation level (h) (e.g. (e.g. ' + stn.split('_')[0] +'_1)')
                                sys.exit(1)

                        elif stn.split('_')[0] == 'autocorr':
                            if not len(stn.split('_'))==3:
                                raise Exception ('The statistic <' + stn.split('_')[0] +  '> has to include the lag (l) and the aggregation level (h) (e.g. autocorr_1_1 (autocorr_l_h))')
                                sys.exit(1)
                            else:
                                if not str(stn.split('_')[1]).isnumeric():
                                    raise Exception ('Aggregation level (h) fot the statistic <' + stn.split('_')[0] +  '> has to be an integer')
                                    sys.exit(1)
                                if not str(stn.split('_')[1]).isnumeric():
                                    raise Exception ('Lag (l) fot the statistic <' + stn.split('_')[0] +  '> has to be an integer')
                                    sys.exit(1)
            else:
                if '_' in stn:
                    raise Exception ('The statistic <mean> should not includ any king of Lag (l) or aggregation level (h). The correct form is <mean>')
                    sys.exit(1)



        if not len(self.statistics_name)==len(self.weights):
            raise Exception ('Hyperparameters <statistics_name> and <weights> have to have the same lenght')
            sys.exit(1) 

        
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
                sys.exit(1)

            res_tuple = type(Seasonality[0]) is tuple
            if res_tuple == False:
                Seasonality_flat_sort = np.sort(Seasonality)
            else:
                Seasonality_flat_sort = np.sort([item for sublist in Seasonality for item in sublist])
            if not len(Seasonality_flat_sort)==12:
                raise Exception ('The longitud of the hyperparameter <Seasonality_user> have to be equal to 12.')
                sys.exit(1)
    
            else:
                if not all(Seasonality_flat_sort == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]):
                    raise Exception ('All the months of the year (from 1 to 12) should be included in the hyperparameter <Seasonality_user>')
                    sys.exit(1)
                
        else:
            raise Exception ('The name of the seasonality type is not correct, check that the name is among those listed below: annual, monthly, seasonal, user_defined')
            sys.exit(1)       
        self.Seasonality_str = Seasonality_str        
        self.Seasonality = Seasonality
        
        
            
        if self.temporal_resolution != 'h' and self.temporal_resolution != 'd':
            raise Exception ("The time resolution of the point model is poorly defined. It should be 'h' or 'd'.")
            sys.exit(1) 
        
        if self.process != 'storms' and self.process != 'normal':
            raise Exception ("The process type is wrongly defined, it should be 'storms' or 'normal'.")
            sys.exit(1) 
            
            
        statistic_list = ['mean_h', 'var_h', 'autocorr_l_h','fih_h', 'fiWW_h', 'fiDD_h', 'M3_h']
        
        for i in self.statistics_name:
            if i.split('_')[0]== 'autocorr':
                st =  i.split('_')[0]+'_'+'l'+'_'+'h'
                if st not in statistic_list:
                    raise Exception ('The '+i+' statistic has not been correctly defined.')
                    sys.exit(1) 
            else:
                st = i.split('_')[0]+'_'+'h'
                if st not in statistic_list:
                    raise Exception ('The '+i+' statistic has not been correctly defined.')
                    sys.exit(1) 
            
            
        
        
        
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

        statistics_name_string = ['mean', 'var', 'autocorr','fih', 'fiWW', 'fiDD', 'M3']
        temporal_resolution_string = ['d', 'h']
        process_string = ['normal', 'storms']

        if self.process not in process_string:
            raise Exception ('The name of the hyperparameter <process> is not correct, check that the name is among those listed below: normal, process')
            sys.exit(1)  

        if self.temporal_resolution not in temporal_resolution_string:
            raise Exception ('The name of the hyperparameter <temporal_resolution> is not correct, check that the name is among those listed below: h, d')
            sys.exit(1) 

        for stn in self.statistics_name:
            if stn.split('_')[0] not in statistics_name_string:
                raise Exception ('The name of the statistic <' + stn.split('_')[0] +  '> is not correct, check that the name is among those listed below: mean_h', 'var_h', 'autocorr_l_h','fih_h', 'fiWW_h', 'fiDD_h', 'M3_h')
                sys.exit(1) 

        if isinstance(self.year_ini, str) or isinstance(self.year_fin, str):
            raise Exception ('The hyperparameter <year_fin> and <year_ini> have to be a integer not a string')
            sys.exit(1) 
        else:
            if str(self.year_ini).isnumeric() and str(self.year_fin).isnumeric():
                if not self.year_fin > self.year_ini:
                    raise Exception ('The hyperparameter <year_fin> has to be higher than the hyperparameter <year_ini>')
                    sys.exit(1) 
            else:
                raise Exception ('The hyperparameter <year_fin> and <year_ini> have to be a integer')
                sys.exit(1) 

        for stn in self.statistics_name:
            if  not 'mean' in stn:
                if not '_' in stn:
                    raise Exception ('Lag (l) or aggregation level (h) is not included for the statistic <' + stn + '>. (e.g. ' + stn +'_1)')
                    sys.exit(1) 
                else:
                    if stn.split('_')[0] not in statistics_name_string:
                        raise Exception ('The name of the statistic <' + stn.split('_')[0] +  '> is not correct, check that the name is among those listed below: mean_h', 'var_h', 'autocorr_l_h','fih_h', 'fiWW_h', 'fiDD_h', 'M3_h')
                        sys.exit(1) 
                    else:
                        if stn.split('_')[0] in ['var','fih', 'fiWW', 'fiDD', 'M3']:
                            if len(stn.split('_'))==2:
                                if not str(stn.split('_')[1]).isnumeric():
                                    raise Exception ('Aggregation level (h) fot the statistic <' + stn.split('_')[0] +  '> has to be a integer (e.g. ' + stn +'_1)')
                                    sys.exit(1)
                            else:
                                raise Exception ('The statistic <' + stn.split('_')[0] +  '> only accept  the aggregation level (h) (e.g. (e.g. ' + stn.split('_')[0] +'_1)')
                                sys.exit(1)

                        elif stn.split('_')[0] == 'autocorr':
                            if not len(stn.split('_'))==3:
                                raise Exception ('The statistic <' + stn.split('_')[0] +  '> has to include the lag (l) and the aggregation level (h) (e.g. autocorr_1_1 (autocorr_l_h))')
                                sys.exit(1)
                            else:
                                if not str(stn.split('_')[1]).isnumeric():
                                    raise Exception ('Aggregation level (h) fot the statistic <' + stn.split('_')[0] +  '> has to be an integer')
                                    sys.exit(1)
                                if not str(stn.split('_')[1]).isnumeric():
                                    raise Exception ('Lag (l) fot the statistic <' + stn.split('_')[0] +  '> has to be an integer')
                                    sys.exit(1)
            else:
                if '_' in stn:
                    raise Exception ('The statistic <mean> should not includ any king of Lag (l) or aggregation level (h). The correct form is <mean>')
                    sys.exit(1)
        
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
                sys.exit(1) 

            res_tuple = type(Seasonality[0]) is tuple
            if res_tuple == False:
                Seasonality_flat_sort = np.sort(Seasonality)
            else:
                Seasonality_flat_sort = np.sort([item for sublist in Seasonality for item in sublist])
            if not len(Seasonality_flat_sort)==12:
                raise Exception ('The longitud of the hyperparameter <Seasonality_user> have to be equal to 12.')
                sys.exit(1)
    
            else:
                if not all(Seasonality_flat_sort == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]):
                    raise Exception ('All the months of the year (from 1 to 12) should be included in the hyperparameter <Seasonality_user>')
                    sys.exit(1)
                
        else:
            raise Exception ('The name of the seasonality type is not correct, check that the name is among those listed below: annual, monthly, seasonal, user_defined')
            sys.exit(1) 
                
        self.Seasonality_str = Seasonality_str        
        self.Seasonality = Seasonality