# NEOPRENE: Neyman-Scott Process Rainfall Emulator

The **NEOPRENE** library implements a rectangular pulses model for rainfall emulation based on the Neyman-Scott process. The emulator may be used to generate synthetic rainfall time series that reproduce observed statistics at different temporal aggregations. It has been designed with rainfall dissaggregation and extreme rainfall analysis in mind.

The description of the Neyman-Scott Process -or Spatio-temporal Neyman-Scott Rectangular Pulses Model- can be found in the [doc folder](doc).

A paper describing the library has been sent for review to _Environmental Modelling & Software_.

Other papers by the authors where -previous incarnations of- the **NEOPRENE** library has been used and the mathematical model has been described are:

+ Diez-Sierra, J.; del Jesus, M. Subdaily Rainfall Estimation through Daily Rainfall Downscaling Using Random Forests in Spain. Water **2019**, _11_, 125. [https://doi.org/10.3390/w11010125](https://doi.org/10.3390/w11010125)
+ del Jesus, M.; Rinaldo, A.; Rodriguez-Iturbe, I. Point rainfall statistics for ecohydrological analyses derived from satellite integrated rainfall measurements. Water Resources Research **2015**, _51(4)_, 2974-2985. [https://doi.org/10.1002/2015WR016935](https://doi.org/10.1002/2015WR016935)

## Contents

| Directory | Contents |
| :-------- | :------- |
|  [NSRP](NSRP) | Python code to calibrate the NSRPM (Neyman-Scott Rectangular Pulse Model) and simulate synthetic rainfall series.
|  [STNSRP](STNSRP) | Python code for calibrate the STNSRPM (Spatio-Temporal Neyman-Scott Rectangular Pulse Model) and simulate multisite rainfall series (in progress).
| [doc](doc) | Description of the model.
| [notebooks](notebooks) |  Jupyter notebooks with examples on how to calibrate, simulate and validate a Neyman-Scott model using the library. Examples on how to perform a daily-to-hourly rainfall disaggregation using the synthetic series is also included.

## Requirements

Scripts and (jupyter) notebooks are provided in [Python](https://www.python.org/) to ensure reproducibility and reusability of the results. The simplest way to match all these requirements is by using a dedicated [conda](https://docs.conda.io) environment, which can be easily installed by issuing:

```sh
conda create -n NEOPRENE pip jupyter
conda activate NEOPRENE
pip install NEOPRENE
```

## Examples of use

Examples of use of the `NEOPRENE` library are available in the form of [jupyter notebooks](notebooks). To run the examples follow the following steps:

1. Download the folder [notebooks](notebooks) from the github repository, or navigate to the folder should you have cloned the repo.
2. Open jupyter notebook of Jupyter Lab (type `jupyter notebook` or `jupyter lab`  in the terminal)
3. Open one of the tests available in the [notebooks](notebooks) folder with jupyter notebook  (e.g. [NSRP_test.ipynb](https://github.com/IHCantabria/NEOPRENE/blob/main/notebooks/NSRP_test.ipynb))

## Errata and problem reporting

To report an issue with the library, please fill a GitHub issue.

## Contributors

The original version of the library was developed by:

+ Javier Díez-Sierra
+ Salvador Navas
+ Manuel del Jesus
