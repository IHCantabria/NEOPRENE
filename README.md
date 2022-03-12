# NEOPRENE: Neyman-Scott Process Rainfall Emulator

[![DOI](https://zenodo.org/badge/409946207.svg)](https://zenodo.org/badge/latestdoi/409946207)

The **NEOPRENE** library implements a rectangular pulses model for rainfall emulation based on the Neyman-Scott process. The emulator may be used to generate multi-site synthetic rainfall time series that reproduce observed statistics at different temporal aggregations. It has been designed with rainfall dissaggregation and extreme rainfall analysis in mind.

The description of the Neyman-Scott Process -or Space-time Neyman-Scott Rectangular Pulses Model- can be found in the [doc folder](https://github.com/IHCantabria/NEOPRENE/tree/main/doc).

A paper describing the library has been sent for review to _Environmental Modelling & Software_.

Other papers by the authors where -previous incarnations of- the **NEOPRENE** library has been used and the mathematical model has been described are:

+ Diez-Sierra, J.; del Jesus, M. Subdaily Rainfall Estimation through Daily Rainfall Downscaling Using Random Forests in Spain. Water **2019**, _11_, 125. [https://doi.org/10.3390/w11010125](https://doi.org/10.3390/w11010125)
+ del Jesus, M.; Rinaldo, A.; Rodriguez-Iturbe, I. Point rainfall statistics for ecohydrological analyses derived from satellite integrated rainfall measurements. Water Resources Research **2015**, _51(4)_, 2974-2985. [https://doi.org/10.1002/2015WR016935](https://doi.org/10.1002/2015WR016935)

## Test the library

If you are curious about how the library works or what it can do, I invite you to go to the **Releases** section of this webpage (on the right-hand side of the page) and download the executable file **NEOPRENE-Setup** for your operative system. This executable file will check if Jupyterlab Desktop is installed in your computer. If it is not, it will download the installation program for you to install Jupyterlab. After Jupyterlab is installed, **NEOPRENE-Setup** will launch the example notebooks for you. Then you can test the library and check its functionality in action.

## Contents

| Directory | Contents |
| :-------- | :------- |
|  [NSRP](https://github.com/IHCantabria/NEOPRENE/tree/main/NEOPRENE/NSRP) | Python code to calibrate the NSRPM (Neyman-Scott Rectangular Pulse Model) and simulate single-site synthetic rainfall series.
|  [STNSRP](https://github.com/IHCantabria/NEOPRENE/tree/main/NEOPRENE/STNSRP) | Python code for calibrate the STNSRPM (Space-Time Neyman-Scott Rectangular Pulse Model) and simulate multi-site synthetic rainfall series.
| [doc](https://github.com/IHCantabria/NEOPRENE/tree/main/doc) | Description of the model.
| [notebooks](https://github.com/IHCantabria/NEOPRENE/tree/main/notebooks) |  Jupyter notebooks with examples on how to calibrate, simulate and validate a Neyman-Scott model using the library. Examples on how to perform a daily-to-hourly rainfall disaggregation using the synthetic series are also included.

## Requirements

Scripts and (jupyter) notebooks are provided in [Python](https://www.python.org/) to ensure reproducibility and reusability of the results. The simplest way to match all these requirements is by using a dedicated [conda](https://docs.conda.io) environment, which can be easily installed by issuing:

```sh
conda create -n NEOPRENE pip jupyter
conda activate NEOPRENE
pip install NEOPRENE
```

## Examples of use

Examples of use of the `NEOPRENE` library are available in the form of [jupyter notebooks](https://github.com/IHCantabria/NEOPRENE/tree/main/notebooks). To run the examples follow the following steps:

1. Download the folder [notebooks](https://github.com/IHCantabria/NEOPRENE/tree/main/notebooks) from the github repository, or navigate to the folder should you have cloned the repo.
2. Open jupyter notebook of Jupyter Lab (type `jupyter notebook` or `jupyter lab`  in the terminal)
3. Open one of the tests available in the [notebooks](https://github.com/IHCantabria/NEOPRENE/tree/main/notebooks) folder with jupyter notebook  (e.g. [NSRP_test.ipynb](https://github.com/IHCantabria/NEOPRENE/blob/main/notebooks/NSRP_test.ipynb), [STNSRP_test.ipynb](https://github.com/IHCantabria/NEOPRENE/blob/main/notebooks/STNSRP_test.ipynb))

## Errata and problem reporting

To report an issue with the library, please fill a GitHub issue.

## Contributors

The original version of the library was developed by:

+ Javier Diez-Sierra
+ Salvador Navas
+ Manuel del Jesus

## License

Copyright 2021 Instituto de Hidráulica Ambiental "IHCantabria". Universidad de Cantabria.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
