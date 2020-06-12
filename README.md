# CT-sim

This project utilizes the package built by Jessica Coyle to investigate the relationship between detection, landscape similarity, and temporal occupancy with simulated species. 

Package repository: https://github.com/hurlbertlab/core-transient-simulation
Full output: XXXX

## Directory structure
The directory structure is nested because multiple experiments were run using this package, so the nestedness of the master repository has been preserved to demonstrate how this experiment (EXP4) fits into the larger core-transient-simulation repository. However, the manuscript only covers the analysis of a single experiment, thus why this repository is distinct from the master repository.

- Code
  - Parameters: nested folder containing the parameters used to set the simulation's initial and summary parameters for each set of dispersal kernels
  - Scripts: visualization scripts for the main analysis (disperal kernel 4), as well as supplementary material (dispersal kernel 2, 8). The script summarizes the output from the simulation in order to create the plots in the manuscript.

- Results
  - EXP4 - examples: A single example output file for one run of the simulation using the conditions disperal kernel 4 and habitat heterogeneity 0.5. Each set of conditions was run for 50 runs. See above for location of the full output.
  - Summary - examples: A single summary file generated from the visualize simulation script, at both the pixel level across species (d-g4_hp-0.5_summary), and for each species (pixel_xclass_summary_bysp). Each set of conditions was summarized, see above for location of the full output.
  - Plots - plots created by the visualization script using the summarized output, in PDF form. 
  
