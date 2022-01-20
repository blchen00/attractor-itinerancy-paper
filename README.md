# attractor-itinerancy-paper

Matlab codes for generating figures in the paper: 

B. Chen, P. Miller, "Attractor-state itinerancy in neural circuits with synaptic depression", *Journal of Mathematical Neuroscience 10, 15 (2020)*.

The repository contains .m files to generate figures in the paper. 

"n1_phase_diagram.m" 
It generates phase diagrams for a single unit's final states with different amplitudes and durations.

"n2_phase_diagram.m" 
It generates phase diagrams for two units' final states with different amplitudes and durations.

"n2_time_evolution"
It computes the evolution of a two-unit system under constant stimuli. The state variables are the firing rate r, the synaptic input s, and the depression d. 

"n2_basins_wij.m"
It estimates the basins of attraction of stable fixed points in a two-unit system, with the cross-coupling weights wij are changing. 

"n5_attractors.m"
It computes the total number of attractors and the number of reached final states in a five-unit system with random cross-couplings. To plot the data, please run "n5_plot.m".
