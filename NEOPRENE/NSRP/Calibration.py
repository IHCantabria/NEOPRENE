from NEOPRENE.NSRP.MathematicalPropertiesNSRP import *
from NEOPRENE.NSRP.utils import *
from NEOPRENE.NSRP.libs_NSRP import *
import time
import yaml
import ast
from NEOPRENE.NSRP.parameters_calibration import parameters_calibration

import cma
import pandas as pd
import numpy as np


class Calibration(object):
    def __init__(self, hiperparams):
        self.hiperparams = hiperparams

    def _initialize_limits(self):
        if self.hiperparams.temporal_resolution == 'd':
            t = 24
        else:
            t = 1

        lim = np.array([
            [(1 / self.hiperparams.time_between_storms[1]) * t, (1 / self.hiperparams.time_between_storms[0]) * t],
            [self.hiperparams.number_storm_cells[0], self.hiperparams.number_storm_cells[1]],
            [(1 / self.hiperparams.cell_duration[1]) * t, (1 / self.hiperparams.cell_duration[0]) * t],
            [(1 / self.hiperparams.cell_intensity[1]) * t, (1 / self.hiperparams.cell_intensity[0]) * t],
            [(1 / self.hiperparams.storm_cell_displacement[1]) * t,
             (1 / self.hiperparams.storm_cell_displacement[0]) * t]])
        return lim

    def _create_dataframes(self, lim):
        if self.hiperparams.process == 'normal':
            param_s = ['landa', 'ipsilon', 'eta', 'xi', 'betha']
            Dataframe_params = pd.DataFrame(index=param_s, columns=self.hiperparams.Seasonality_str)
            n_p = 5
        elif self.hiperparams.process == 'storms':
            param_s = ['landa1', 'ipsilon1', 'eta1', 'xi1', 'betha1', 'landa2', 'ipsilon2', 'eta2', 'xi2', 'betha2']
            Dataframe_params = pd.DataFrame(index=param_s, columns=self.hiperparams.Seasonality_str)
            n_p = 10

        limites = np.array(lim)
        return limites, Dataframe_params, n_p

    def _calibrate_parameters(self, limites, Dataframe_params, n_p):
        Error_dataframe = pd.DataFrame(index=self.hiperparams.Seasonality_str)
        Error_dataframe['Total Error'] = np.nan

        self.statististics_fit_dataframe = pd.DataFrame(index=self.hiperparams.Seasonality_str)

        for pri, prii in enumerate(self.hiperparams.Seasonality):
            evaluator = evaluateInd_PSO(self.statististics_dataframe.loc[:, str(prii)].values, self.hiperparams.weights,
                                        self.hiperparams.process, self.hiperparams.statistics_name)

            def objective_function(params):
                return evaluator.totalError(params)

            # Inicializar CMA-ES
            sigma = 0.5
            initial_guess = [(low + high) / 2 for low, high in limites]
            bounds = [limites[:, 0], limites[:, 1]]  # Formato adecuado para CMA-ES

            es = cma.CMAEvolutionStrategy(initial_guess, sigma, {'bounds': bounds})

            # Ejecutar optimización
            es.optimize(objective_function)

            # Obtener resultados
            Fitted_parameters = es.result.xbest
            total_error = es.result.fbest

            # Almacenar resultados en el DataFrame
            sta_fit = evaluator.compute(Fitted_parameters)
            self.statististics_fit_dataframe.loc[:, str(prii)] = sta_fit[0]
            Error_dataframe.loc[str(prii)] = total_error

            # Actualizar DataFrame de parámetros ajustados
            if self.hiperparams.process == 'normal':
                param_v = Fitted_parameters[:5]
            elif self.hiperparams.process == 'storms':
                param_v = Fitted_parameters

            Dataframe_params.loc[:, str(prii)] = param_v

        resuls = parameters_calibration(Dataframe_params, self.statististics_fit_dataframe,
                                        self.statististics_dataframe, Error_dataframe)
        return resuls

    def __call__(self, statistics, verbose=False):
        self.statististics_dataframe = statistics.statististics_dataframe
        lim = self._initialize_limits()
        limites, Dataframe_params, n_p = self._create_dataframes(lim)
        results = self._calibrate_parameters(limites, Dataframe_params, n_p)
        return results
    
