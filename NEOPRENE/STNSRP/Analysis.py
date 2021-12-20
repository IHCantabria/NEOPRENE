'''
Library containing classes for calibrating model parameters.

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
    
'''

from NEOPRENE.STNSRP.utils import *

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
        
        
    

