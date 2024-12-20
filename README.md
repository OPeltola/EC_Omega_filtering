# EC_Omega_filtering
These codes can be used to reproduce the results presented in Peltola et al. (2025): "Towards an enhanced metric for detecting vertical flow decoupling in eddy covariance flux observations" published in Agricultural and Forest Meteorology.
## Installation
The codes have been run under Windows OS with Python 3.11.6. Install the required Python packages by running 
```Shell
pip install -r requirements.txt
```
## Usage
### Reproducing figures
You can reproduce the figures shown in the manuscript by downloading the data from Zenodo repository (<https://doi.org/10.5281/zenodo.14259134>) and then by running
```Shell
python plot_figures.py
```
Note that the data needs to be located in subfolder `data/`.
### Processing data
The data files included in the Zenodo repository can be recreated by downloading data from NEON data portal and ICOS carbon portal, see DOI for each dataset in the manuscript, and then by running
```Shell
python process_data.py
```
Note that `.env` file is needed, there the location of each dataset on the computer should be described.
### Example Omega calculation
For calculating Omega decoupling parameter for a few example cases, you can run
```Shell
python example_omega_calculation.py
```


## Citing
If you are using this codebase, please cite: Peltola, O., Aslan, T., Aurela, M., Lohila, A., Mammarella, I., Papale, D., Thomas, C.K., Vesala, T., Laurila, T., 2025. Towards an enhanced metric for detecting vertical flow decoupling in eddy covariance flux observations. Agricultural and Forest Meteorology 362, 110326. https://doi.org/10.1016/j.agrformet.2024.110326

