# C-fiber Model with Waiting-State Sodium Channel Inactivation

This computational model simulates the biophysical properties of C-fibers, a class of unmyelinated sensory nerve fibers involved in pain perception. The model incorporates detailed ion channel dynamics, including sodium and potassium channels.

In this branch, the fast inactivation dynamics of the sodium channels Nav1.7, Nav1.8, and Nav1.3 are modified using a multi-state waiting-state model based on Köster et al. (2025):

h0 -> h1 -> hw1 -> hw2 -> hw3 -> hw4 -> hw5 -> h0

The effective channel availability is calculated as:

h_eff = 1 - h0

For each modified sodium channel, the current-producing mechanism and the
fast-inactivation mechanism are implemented separately. The inactivation
state is passed to the corresponding current mechanism using a NEURON POINTER.

## Requirements:
- NEURON v7.8 or higher
- Python
- NumPy, pandas, and Matplotlib for data analysis and plotting

## Main Modifications
### Model Files:
- `main.py`: runs the C-fiber model and records membrane potential, sodium currents, and channel-state variables
- `defineCell.py`: inserts the additional inactivation mechanisms and connects
  them to the corresponding sodium-channel mechanisms

### Sodium Channel mechanisms
- `nattxs.mod` and `nattxs_h.mod`: modified Nav1.7 fast inactivation
- `DNav18.mod` and `DNav18_h.mod`: modified Nav1.8 fast inactivation
- `nav13.mod` and `nav13_h.mod`: modified Nav1.3 fast inactivation

### Plotting:
- `plot.py`: compares ADS and recovery-cycle results between the original and modified models
- `plot_cycle.py`: compares recovery cycles at different baseline frequencies

## Citation:
If you use this model in your research, please cite the original model publication:

"Maxion, A., et al. (2023). A modeling study to dissect the potential role of voltage-gated ion channels in activity-dependent conduction velocity changes as identified in small fiber neuropathy patients. Frontiers in Computational Neuroscience, 17. https://doi.org/10.3389/fncom.2023.1265958"

"Köster, Phil Alexander, et al. (2025). Nociceptor sodium channels shape subthreshold phase, upstroke, and shoulder of action potentials. Journal of General Physiology 157.2 : e202313526. https://doi.org/10.1085/jgp.202313526"

## License:
This project is licensed under the Apache License 2.0.
