## This is the GitHub Repository for the publication: 

# Biomass-derived Diformylxylose as a Renewable Solvent for Biocatalysis Applications

In this repository we provide the code that was used for MD simulations in diformylxylose (DFX) using OpenMM [1] and analysis of enzyme stability related parameters using MDTraj [2] and EvoEF2 [3].

[1] Eastman, P., Galvelis, R., Peláez, R. P., Abreu, C. R., Farr, S. E., Gallicchio, E., ... & Markland, T. E. (2023). OpenMM 8: molecular dynamics simulation with machine learning potentials. The Journal of Physical Chemistry B, 128(1), 109-116.

[2] McGibbon, R. T., Beauchamp, K. A., Harrigan, M. P., Klein, C., Swails, J. M., Hernández, C. X., ... & Pande, V. S. (2015). MDTraj: a modern open library for the analysis of molecular dynamics trajectories. Biophysical journal, 109(8), 1528-1532.

[3] Huang, X., Pearce, R., & Zhang, Y. (2020). EvoEF2: accurate and fast energy function for computational protein design. Bioinformatics, 36(4), 1135-1142.


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
Navigate to KRED_36 folder
```bash
cd MD_sims/KRED_36
```
Run the MD simulation of KRED_36
```bash
python MD_sims/KRED_36/MD_simulation.py
```
The resulting trajectories will be found in subfolder trajectories_KRED_36


Navigate to TeSADH_W110A folder
```bash
cd ../TeSADH_W110A
```
Run the MD simulation of TeSADH_W110A
```bash
python MD_simulation.py
```
The resulting trajectories will be found in subfolder trajectories_TeSADH_W110A

## How to analyze the RoG of the trajectories
Activate conda environment
```bash
conda activate dfx_sim
```
Run analysis for KRED_36
```bash
python RoG_analysis/RoG_analysis.py --input_folder MD_sims/KRED_36/trajectories_KRED_36 --output_folder RoG_analysis/KRED_36
```
Run analysis for TeSADH_W110A
```bash
python RoG_analysis/RoG_analysis.py --input_folder MD_sims/TeSADH_W110A/trajectories_TeSADH_W110A --output_folder RoG_analysis/TeSADH_W110A
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
python EvoEF2_batch.py --input_folder EvoEF2_inputs --output_folder EvoEF2_outputs
 
```

# References

If you utilize this code, please cite (as soon as published):
Fatma Feyza Özgen et al. Biomass-derived Diformylxylose as a Renewable Solvent for Biocatalysis Applications (submitted)
