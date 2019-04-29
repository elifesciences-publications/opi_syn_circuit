# opi_syn_circuit

Analysis and data repository for the opi_syn_circuit project:
Related to Birdsong*, Jongbloets*, Engeln, Wang, Scherrer, Mao, 2019, eLife, Synapse-specific opioid modulation of thalamo-cortico-striatal circuits. *These authors contributed equally*

First look through figure_data_script_overview.csv. This file links all the figure panels in the main text and rebuttal to the source data and the scripts in R and Matlab.

Without the need for installing R or Matlab you can open either
'opi_syn_circuit_plots_R.html' or 'opi_syn_circuit_plots_R.pdf' to browse through the panels and statistics. If you want to actually go through the code and do meta-analysis, please look at the 'opi_syn_circuit_plots_R.ipynb'. This file contains the jupyter notebook source for the html and pdf files.

Please be aware that for the R scripts certain packages are required. Therefore read the jupyter notebook or R scripts package import commands to find the required packages.

Folders:
data: contains source data and preprocessed data for analysis, plotting,        statistics
importaxo: contains matlab script to read axograph files into matlab
movingmean: contains matlab script to smooth waveform traces using a moving average.
