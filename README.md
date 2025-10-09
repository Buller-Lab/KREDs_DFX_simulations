## This is the GitHub Repository for the publication: 

# Biomass-derived Diformylxylose as a Renewable Solvent for Biocatalysis Applications

In this repository we provide the code that was used for MD simulations in diformylxylose (DFX) using OpenMM [1] and analysis of enzyme stability related parameters using MDTraj and EvoEF2.[2,3] 

[1] Guerra, João Victor da Silva, et al. "pyKVFinder: an efficient and integrable Python package for biomolecular cavity detection and characterization in data science." BMC bioinformatics 22.1 (2021): 607.

[2] Wohlwend, Jeremy, et al. "Boltz-1 democratizing biomolecular interaction modeling." BioRxiv (2025): 2024-11.

[3] Passaro, Saro, et al. "Boltz-2: Towards accurate and efficient binding affinity prediction." BioRxiv (2025): 2025-06.

[4] McGibbon, Robert T., et al. "MDTraj: a modern open library for the analysis of molecular dynamics trajectories." Biophysical journal 109.8 (2015): 1528-1532.

[5] Pedregosa, Fabian, et al. "Scikit-learn: Machine learning in Python." the Journal of machine Learning research 12 (2011): 2825-2830.

# Installation

We recommend to run this code on UNIX based systems such as Ubuntu. This repository can be downloaded to your local machine via the command:
```bash
git clone https://github.com/Buller-Lab/KREDs_DFX_simulations
```
then navigate into the cloned repository with:
```bash
cd KREDs_DFX_simulations
```

this should only take a few seconds.

# System Requirements

## Hardware requirements

This code was developed and tested on the following hardware:

- CPU: AMD Ryzen Threadripper 3970X 32-Core Processor
- Memory: 130 GiB RAM
- GPU: 2x NVIDIA GeForce RTX 3090

## Software requirements
To create conda environments with necessary dependencies, run:
```bash
conda env create --file dfx_sim.yml
```
# Instructions for use
## The followings scripts are provided:
- MD_simulation.py (Custom MD simulation script for KRED solvent simulations)
- RoG_analysis.py (Radius of Gyration analysis of all trajectories)
- EvoEF2_batch.py (Total energy calculation of final simulation frames) 
## How to run the MD simulations (recommended in screen or on HPC)
Activate conda environment
```bash
conda activate dfx_sim
```
Run the MD simulation of KRED_36
```bash
python MD_sims/KRED_36/MD_simulation.py
```
The resulting trajectories will be found in subfolder trajectories_KRED_36

Run the MD simulation of TaSADH_W110A
```bash
python MD_sims/TaSADH_W110A/MD_simulation.py
```
The resulting trajectories will be found in subfolder trajectories_TaSADH_W110A

## How to analyze the RoG of the trajectories
Activate conda environment
```bash
conda activate dfx_sim
```

## How to calculate the EvoEF2 total energy 
Activate conda environment
```bash
conda activate dfx_sim
```
Navigate in working directory
```bash
 cd ddG_calc
```
Clone EvoEF2 repository
```bash
git clone https://github.com/tommyhuangthu/EvoEF2.git
```
Navigate to EvoEF2 folder
```bash
cd EvoEF2
```
Make EvoEF2 build executable
```bash
chmod +x build.sh
```
Install EvoEF2
```bash
bash build.sh
```
Navigate back to ddG_calc folder
```bash
cd ..
```
Run the energy calculation if the input folder in one batch generating results in the output folder
```bash
python EvoEF2_batch.py --input_folder EvoEF2_inputs --output_colder EvoEF2_outputs
 
```

# References

If you utilize this code, please cite:
Daniela Schaub, Alice Lessing et al. A tailored enzyme cascade facilitates DNA-encoded library technology and gives access to a broad substrate scope, 23 September 2025, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-7598475/v1]
