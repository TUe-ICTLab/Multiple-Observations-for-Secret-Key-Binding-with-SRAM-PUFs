# SRAMPUF-MultipleObservations
MATLAB code used for generating the plots in our "Multiple Observations for Secret-Key Binding with SRAM PUFs" paper.

L. Kusters and F.M.J. Willems, "Multiple Observations for Secret-Key Binding with SRAM PUFs," in preprints, doi: https://doi.org/10.20944/preprints202104.0229.v1

Files description:

- calc_mutinf.m : generate Fig. 2 and 4 from the paper
- example_simulations_call : example of calling the simulations
- example_plot_results : example of plotting simulation results for Fig. 6 and 9
- /sim_results : example results for running the simulations
- /functions/sim_MOhds.m : Monte Carlo simulation of Multiple Observations scheme
- /functions/sim_SMOhds.m: Monte Carlo simulation of Sequential Multiple Observations scheme


This package requires the following MATLAB toolboxes
- Communications toolbox
- 5G Toolbox
- DSP System Toolbox
- Statistics and Machine Learning Toolbox
