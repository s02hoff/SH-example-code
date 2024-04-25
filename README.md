# SH-example-code

#below are descriptions for each analysis file 

1) IslandMYSE_WNS: here we analyze disease dynamics for populations of northern myotis, including pathogen prevalence and infection intensity between island and mainland populations (using GLMMs), population growth rate data (to measure declines at hibernating colonies), comparison of spring bat weights in pre- and post-white-nose syndrome years (LM), prevalence and contamination (load) of the environmental reservoir (GLM), and capture probabilities at both locations to model the decline in probability of capturing a bat on a single night (binomial GLM)

2) create_capt_hist: this file loads in our presence/absence data and environmental covariates, transforms data into capture history, normalizes covariates, and creates a capture history imp file that can be read by program Mark
   
3) run_bat_mods: this file converts the data and creates the model structure, and runs various multi-scale single season occupancy models to assess covariate relationships with landscape occupancy (psi), local occupancy (theta), and detection (p)

4) oSCR_workflow: here we conduct a spatial capture-recapture analysis to determine the effect of environmental disturbance on species abundance and density
   
5) KF_analysis: here we conduct known-fate analysis from telemetry data to estimate survival over time and by individual parameters (age, weight, etc.)
   
6) amakini_density: This analysis estimates density of animals from data collected with a distance sampling method
