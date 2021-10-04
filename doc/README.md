# The NEOPRENE library

## Concept

The **NEOPRENE** library allows to work with Neyman-Scott models for rainfall emulations in an easy way. The library has two main functionalities:

1. **Model Calibration:** Taking daily or hourly time series, or their characteristics statistics, a set of optimal parameters are computed. These parameters allow to reproduce the main statistical properties of the time series, or at least, the provided statistics.

   ![Scheme of the calibration process](Calibration.png)

2. **Rainfall Generation:** Given a set of parameters, time series at the daily or hourly scale, or their characteristics statistics can be generated.

   ![Scheme of the generation process](Generation.png)

A mathematical description of the model is given in the [linked document](MathematicalDescription.md).

The following sections present the calibraton and validation hyperparameters, as well as some elements that should be considered when using the **NEOPRENE** library.

## Implementation

Currently, the library can only represent rainfall at a given point. The basic model used is represented by the following parameters:

1. Storm arrivals follow a **Poisson** process of parameter &lambda; expressed in storms per square kilometer and per day.
2. Each storm is characterized by a number of cells that follow a **Poisson** process of parameter &upsilon; expressed in storm cells per storm.
3. Each storm cell is characterized by three variables:
   1. The cell lag from the storm origin that follows an **Exponential** process of parameter &beta; measured in hours^-1.
   2. The cell duration that follows an **Exponential** process of parameter &epsilon; measured in hours^-1
   3. The cell rainfall intensity that follow an **Exponential** process of parameter &chi; measured in hours per milimeter.

## Calibration hyperparameters

+ **Data:** Pandas ```DataFrame``` that contains the original time series that is to be emulated using **NEOPRENE**.

+ **Seasonality:** Python ```list``` that configures the desired seasonality for the model. Calibration can be done in a monthly basis, by season or by year.
  + _Anual calibraton_: The library assumes that a single set of parameters is able to capture the dynamics for the whole year.

    ```python
    Seasonality=list()
    Seasonality.append((1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ,11, 12))
    ```

  + _Seasonal calibraton_: The library merges the different months into the prescribed groups, fitting a different set of parameters per group.
  
    ```python
    Seasonality=list()
    Seasonality.append((1, 2, 3))
    Seasonality.append((4, 5, 6))
    Seasonality.append((7, 8, 9))
    Seasonality.append((10, 11, 12))
    ```

  + _Monthly calibration_: The library fits a different set of parameters for each month of the time series.
  
    ```python
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
    ```

+ **temporal_resolution:** ```string``` specifying the temporal resolution of the time series provided to the calibration process.

    ```python
    temporal_resolution = 'h' # Hourly resolution
    temporal_resolution = 'd' # Daily resolution
    ```

+ **process**: ```string``` configuring the model type.

  ```python
  process = 'normal' # Only one type of storm is considered
  process = 'storms' # Convective and Frontal storm are considered
  ```

+ **Intensity_function:** ```string``` that configures the distribution used for the rainfall intensity. Currently, only exponential distributions are considered.

  ```python  
  Intensity_function = 'E'
  ```

+ **statistics:** ```list``` of ```strings``` that contain the statistics that have to be considered during the fitting process. The statistics included are:

   ```python
   statistics = ['mean', # Rainfall average
                'var_h', # Variance
                'autocorr_l_h', # Autocorrelation
                'fih_h', # Probability of no rainfall
                'fiWW_h', # Transition probability from rainy period to rainy period
                'fiDD_h', # Transition probability from dry period to dry period
                'M3_h'] # Skewness
   ```

   Statistics may refer to different lags (```l```) and aggregation levels (```h```), where the aggregation levels indicate the number of temporal resolutions over which the value is aggregated.

+ **weights:** ```list``` that contains the weights for computing the total error during the calibration process

  ```python
  weights = [100, 1, 1, 1, 1, 1, 1]
  ```

+ **number_iterations:** ```integer``` to define the maximum number of iterations of the calibration process

    ```python
    number_iterations = 20
    ```

+ **number_bees:** ```integer``` defining the number of particles to be used by the calibration algorithm. Calibration is carried out by means of a Particle Swarm Optimization (PSO) algorithm.

  ```python
  number_bees = 10000
  ```
  
+ **number_initializations:**```integer``` defining the number of initializations to be performed during the calibration procedure.

   ```python
  number_initializations = 1
  ```

+ **time_between_storms:** ```list``` defining the range of inter-storm arrival times. The inter-storm arrival time is the inverse of the &lambda; parameter.

  ```python
  time_between_storms = [20, 1000] # hours
  ```

+ **number_storm_cells:** ```list``` defining the range of possible values for the number of storm cells per storm (&upsilon;)

  ```python
  number_storm_cells = [1, 100] # number of cells per storm
  ```

+ **cell_duration:** ```list``` defining the range of possible storm duration values. It is the inverse of the &epsilon; parameter.

  ```python
  cell_duration = [0.15, 5] # hours
  ```

+ **cell_intensity:** ```list``` defining the range of intensities of each storm cell. It is the inverse of the &chi; parameter.

  ```python
  cell_intensity = [0.5, 2000] # mm / hour 
  ```

+ **storm_cell_displacement:** ```list``` defining the range of acceptable lags between the storm origin and each storm cell. It is the inverse of the &beta; parameter.

  ```python
  storm_cell_displacement = [1.01, 50] # hours
  ```

## Simulation hyperparameters

+ **year_ini and year_fin:** Initial and final year of the simulated time series

```python
year_ini = 2000
year_fin = 2100
```

+ **parameters_simulation**: ```list``` containing the values of the simulated parameters

  ```python
  parameters_simulation = [lanbda, upsilon, beta, chi, epsilon]
  ```

+ **statistics_external**: ```list``` providing the statistics to be reproduced by the generated time series

  ```python
  statistics_external = ['mean', 'var_1', 'var_2', 'var_3', 'var_4', 'autocorr_1_1', 'autocorr_2_1', 'autocorr_3_1','fih_1', 'fiWW_1', 'fiDD_1', 'M3_1']
    ```