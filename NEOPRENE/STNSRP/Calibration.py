# Calibration.py

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
import gc

import numpy as np
import pandas as pd
import logging

np.random.seed(42)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class Calibration(object):
    def __init__(self, hiperparams):
        self.hiperparams = hiperparams

    def __call__(self, statistics, time_series, verbose=False):
        self.statistics_dataframe = statistics.statistics_dataframe.copy()
        self.crosscorr_dataframe = statistics.crosscorr_dataframe.copy()
        self.crosscorr_dataframe_dist = statistics.crosscorr_dist_dataframe.copy()

        t = 24 if self.hiperparams.temporal_resolution == 'd' else 1
        lim = self.calculate_limits(t)

        sta_len = len(self.hiperparams.statistics_name)
        nn_cross = np.sum(['cross' in i for i in self.hiperparams.statistics_name])
        cross_corr_division = 10
        statistics_dataframe_fit = pd.DataFrame(index=self.hiperparams.statistics_name[0:sta_len - nn_cross], columns=self.hiperparams.Seasonality_str)
        crosscorr_dataframe_fit = {key: {key_season: pd.DataFrame(index=self.crosscorr_dataframe[key][key_season].index, columns=self.crosscorr_dataframe[key][key_season].columns)
                                         for key_season in self.crosscorr_dataframe[key].keys()} for key in self.crosscorr_dataframe.keys()}

        param_s, Dataframe_params, n_p, limites = self.setup_parameters(lim)

        logging.info('\n\n' + '#' * 80)
        logging.info('Adjustment of parameters using the Particle Swarm Optimization (PSO)\n\n')

        Error_dataframe = pd.DataFrame(index=self.hiperparams.Seasonality_str)
        Error_dataframe['Total Error'] = np.nan
        statistics_observed_dataframe = pd.DataFrame(index=self.hiperparams.statistics_name)

        for pri, prii in enumerate(self.hiperparams.Seasonality):
            try:
                self.process_seasonality(prii, sta_len, nn_cross, cross_corr_division, statistics_dataframe_fit, crosscorr_dataframe_fit, Error_dataframe, Dataframe_params, n_p, limites, verbose)
            except Exception as e:
                logging.error(f"Failed to process seasonality {prii}: {e}")
                continue

        Dataframe_xi_months = XI_MONTHS(time_series.astype(float), Dataframe_params.astype(float), self.hiperparams.process)
        resuls = parameters_calibration(Dataframe_params, statistics_dataframe_fit, crosscorr_dataframe_fit, self.statistics_dataframe, self.crosscorr_dataframe, self.crosscorr_dataframe_dist, Error_dataframe, Dataframe_xi_months)
        return resuls

    def calculate_limits(self, t):
        if self.hiperparams.storm_radius:
            return np.array([
                [(1 / self.hiperparams.time_between_storms[1]) * t, (1 / self.hiperparams.time_between_storms[0]) * t],
                [self.hiperparams.number_storm_cells[0], self.hiperparams.number_storm_cells[1]],
                [(1 / self.hiperparams.cell_duration[1]) * t, (1 / self.hiperparams.cell_duration[0]) * t],
                [(1 / self.hiperparams.storm_cell_displacement[1]) * t, (1 / self.hiperparams.storm_cell_displacement[0]) * t],
                [(1 / self.hiperparams.cell_radius[1]), (1 / self.hiperparams.cell_radius[0])],
                [(1 / self.hiperparams.storm_radius_p[1]), (1 / self.hiperparams.storm_radius_p[0])]
            ])
        else:
            return np.array([
                [(1 / self.hiperparams.time_between_storms[1]) * t, (1 / self.hiperparams.time_between_storms[0]) * t],
                [self.hiperparams.number_storm_cells[0], self.hiperparams.number_storm_cells[1]],
                [(1 / self.hiperparams.cell_duration[1]) * t, (1 / self.hiperparams.cell_duration[0]) * t],
                [(1 / self.hiperparams.storm_cell_displacement[1]) * t, (1 / self.hiperparams.storm_cell_displacement[0]) * t],
                [(1 / self.hiperparams.cell_radius[1]), (1 / self.hiperparams.cell_radius[0])]
            ])

    def setup_parameters(self, lim):
        if self.hiperparams.process == 'normal':
            if not self.hiperparams.storm_radius:
                param_s = ['landa', 'ipsilon', 'eta', 'betha', 'fi_may']
                Dataframe_params = pd.DataFrame(index=param_s, columns=self.hiperparams.Seasonality_str)
                n_p = 5
                limites = lim[:5]
            else:
                param_s = ['landa', 'ipsilon', 'eta', 'betha', 'fi_may', 'fi_may_s']
                Dataframe_params = pd.DataFrame(index=param_s, columns=self.hiperparams.Seasonality_str)
                n_p = 6
                limites = lim
        elif self.hiperparams.process == 'storms':
            if not self.hiperparams.storm_radius:
                param_s = ['landa1', 'ipsilon1', 'eta1', 'betha1', 'fi_may1', 'landa2', 'ipsilon2', 'eta2', 'betha2', 'fi_may2']
                Dataframe_params = pd.DataFrame(index=param_s, columns=self.hiperparams.Seasonality_str)
                n_p = 10
                limites = np.vstack((lim[:5], lim[:5]))
            else:
                param_s = ['landa1', 'ipsilon1', 'eta1', 'betha1', 'fi_may1', 'fi_may_s1', 'landa2', 'ipsilon2', 'eta2', 'betha2', 'fi_may2', 'fi_may_s2']
                Dataframe_params = pd.DataFrame(index=param_s, columns=self.hiperparams.Seasonality_str)
                n_p = 12
                limites = np.vstack((lim, lim))
        return param_s, Dataframe_params, n_p, limites

    def get_fitted_parameters(self, Fitted_parameters):
        if self.hiperparams.process == 'normal' and not self.hiperparams.storm_radius:
            return Fitted_parameters[:5]
        elif self.hiperparams.process == 'normal' and self.hiperparams.storm_radius:
            return Fitted_parameters[:6]
        elif self.hiperparams.process == 'storms' and not self.hiperparams.storm_radius:
            return Fitted_parameters[:10]
        elif self.hiperparams.process == 'storms' and self.hiperparams.storm_radius:
            return Fitted_parameters[:12]

    def process_seasonality(self, prii, sta_len, nn_cross, cross_corr_division, statistics_dataframe_fit, crosscorr_dataframe_fit, Error_dataframe, Dataframe_params, n_p, limites, verbose):
        cross_corr_stationality_month = {k: self.crosscorr_dataframe[k][prii] for k in self.crosscorr_dataframe.keys()}
        statistics_observed_dataframe = self.statistics_dataframe[prii].mean(axis=1)
        evaluator = evaluateInd_PSO(statistics_observed_dataframe, self.hiperparams.weights, self.hiperparams.process, self.hiperparams.statistics_name,
                                    self.hiperparams.storm_radius, cross_corr_stationality_month)
        if verbose:
            logging.info(f'\nFitting the months = {prii}\n')
        error = 0.001
        error_initializations = []
        parameters_initializations = []
        n_i = 0
        total_error = 10000

        while total_error > error and self.hiperparams.number_initializations > n_i:
            if verbose:
                logging.info(f'Initialization number  = {n_i}')
            n_i += 1
            Steps = 0
            mySwarm = pso.Swarm(n_p, limites, self.hiperparams.number_bees, evaluator, 'min', verbose)

            while True:
                try:
                    if verbose:
                        logging.info(f'Iteration number  = {Steps}')
                    mySwarm.step()
                    Steps += 1
                    b, bp = mySwarm.getBest()
                    current_error = evaluator.totalError(bp)

                    if verbose:
                        logging.info(f'Total error = {current_error}')
                        if total_error <= error:
                            logging.info(f'Stopped calibration since total error is equal or lower than error bound ({error})')

                    if current_error < error or Steps > self.hiperparams.number_iterations:
                        break
                except Exception as e:
                    logging.error(f"Error during optimization: {e}")
                    break

            error_initializations.append(current_error)
            parameters_initializations.append(bp)
            sta_fit = evaluator.compute(bp)
            if len(sta_fit) > 0:  # Check to avoid IndexError
                statistics_dataframe_fit.loc[:, str(prii)] = sta_fit[0:sta_len - nn_cross]
                for nn_ccc, st_cross in enumerate(self.hiperparams.statistics_name[sta_len - nn_cross:]):
                    crosscorr_dataframe_fit[st_cross][prii]['cross_corr'] = sta_fit[(sta_len - nn_cross) + (nn_ccc * cross_corr_division):(sta_len - nn_cross) + ((nn_ccc + 1) * cross_corr_division)]
                    crosscorr_dataframe_fit[st_cross][prii]['dist'] = self.crosscorr_dataframe[st_cross][prii]['dist']
            Fitted_parameters = mySwarm.getBest()[1]
            Error_dataframe.loc[str(prii)] = current_error
            total_error = current_error

        if n_i > 1:
            Error_dataframe.loc[str(prii)] = np.min(error_initializations)
            sta_fit = evaluator.compute(parameters_initializations[np.argmin(error_initializations)])
            if len(sta_fit) > 0:  # Check to avoid IndexError
                statistics_dataframe_fit.loc[:, str(prii)] = sta_fit[0:sta_len - nn_cross]
                for nn_ccc, st_cross in enumerate(self.hiperparams.statistics_name[sta_len:]):
                    crosscorr_dataframe_fit[st_cross][prii]['cross_corr'] = sta_fit[(sta_len - nn_cross) + (nn_ccc * cross_corr_division):(sta_len - nn_cross) + ((nn_ccc + 1) * cross_corr_division)]
                    crosscorr_dataframe_fit[st_cross][prii]['dist'] = self.crosscorr_dataframe[st_cross][prii]['dist']
            Fitted_parameters = parameters_initializations[np.argmin(error_initializations)]

        param_v = self.get_fitted_parameters(Fitted_parameters)
        Dataframe_params.loc[:, str(prii)] = param_v
        if verbose:
            logging.info(f'Best params {param_v}')

        del param_v, Fitted_parameters, evaluator, mySwarm, error_initializations, parameters_initializations, sta_fit
        gc.collect()





    
