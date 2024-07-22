# Phylodynamics Beyond Neutrality

Simulation and analysis code for exploring the effects of purifying selection on viral phylogenies and phylodynamic inference.

For more information, please see the our manuscript:

Koelle, K. and Rasmussen, D.A. Phylodynamics beyond neutrality: The impact of incomplete purifying selection on viral phylogenies and inference. In prep.

## Overview of code

Our simulation study on how purifying selection impacts tree summary statistics can be replicated in full by running **run_slim_sims.py**. This script will run simulations at all unique model paramterizations we consider including each mutational fitness effect s_d, mutation rate and sampling interval considered in the paper.

The **slim** folder contains the scripts used to simulate epidemic trajectories and phylogenetic trees in SLiM.
- **linear_birth_death_fixed_sd.slim** can be run directly in SLiM to simulate trees with a fixed mutational fitness costs s_d.
- **linear_birth_death_low_mu.slim** and **linear_birth_death_high_mu.slim** are intended to be called from within **run_slim_sims.py** as they expect model parameters to be passed from the command line. The mutation rates are fixed and are set at their high/low values as given in the filename.

Note: environment.yml can be used to create a conda environment with SLiM, tskit and other required packages installed:
```
$ conda env create -f environment.yml
```


