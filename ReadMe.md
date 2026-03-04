# Stochastic Computational Model of a C-fiber

This computational model simulates the biophysical properties of C-fibers, a class of unmyelinated sensory nerve fibers 
involved in pain perception. The model incorporates detailed ion channel dynamics, including sodium and potassium channels.
It was translated to Python and adapted to run on a high-performance computing cluster (HPC). 
Stochasticity has been added to the conductances and  gating variables of the three sodium channels: Nav1.7, Nav1.8, and Nav1.9. 

## Requirements:
- NEURON v7.8 or higher
- (For PSO-based optimization) Access to a High-Performance Cluster

## Installation and Usage:
1. Install NEURON: 
  Follow the official installation guide: https://www.neuron.yale.edu/neuron/
2. Compile mod-files:
  Navigate to the MOD_Tigerholm folder and ensure it has all .mod files
  Compile nrnivmodl of MOD_Tigerholm in NEURON
3. Run the Model
   
i) Deterministic   
  To execute the model, use: 
  
  import main  
  main.run()

ii) Stochastic   
  To execute the stochastic model, use: 
  
  import main_noise_Nav  
  main_noise_Nav.run()

iii) Stochastic with mutation implemented   
  To execute the stochastic model with mutation implemented, use: 
  
  import main_noise_Nav_p  
  main_noise_Nav_p.run()

## File Structure:
### Model Files:
- main.py: Creates the nerve cell, runs the simulation, and saves results
- main_noise_Nav.py: Creates the nerve cell, runs the simulation, and saves results (stochastic)
- main_noise_Nav_p.py: Creates the nerve cell, runs the simulation, and saves results (stochasticity and mutation implemented)
- run.py: Example script to run the model
- defineCell.py: Functions for creating the cell 
- defineCell_noise_Nav.py: Functions for creating the cell (stochastic)
- defineCell_noise_Nav_p.py: Functions for creating the cell (stochasticity and mutation implemented)
- dataProcessing.py: Functions for data processing: getFilename(), getData() and calculateLatency()
- stimulationProtocols.py: Predefined stimulation protocols

### Optimization Files:
- particleSwarm2.py, particleSwarm2GEPSO.py, particleSwarm2RPSO.py: Different implementations of the 
Particle Swarm Optimization (PSO) algorithm
- evaluate.py: Functions for evaluating PSO-generated parameter sets

### Plotting:
- plot.py: Functions for visualizing results

### Mechanisms: 
- All ion channel mechanism files (.mod files) are located in MOD_Tigerholm/

## Citation:
If you use this model in your research, please cite the following publication(s):

"Maxion, A., et al. (2023). A modeling study to dissect the potential role of voltage-gated ion channels in 
activity-dependent conduction velocity changes as identified in small fiber neuropathy patients. Frontiers 
in Computational Neuroscience, 17. https://doi.org/10.3389/fncom.2023.1265958"

Roseira, J. 

## License:
This project is licensed under the Apache License 2.0. 
