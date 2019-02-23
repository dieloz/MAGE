# MAGE
MAGE: Multiscale Adaptive Gabor Expansion

This is the README for MATLAB code and scripts associated with the paper "Multiscale Adaptive Gabor Expansion (MAGE): Improved Detection of Transient Oscillatory Burst Amplitude and Phase", by Ryan T. Canolty and Thilo Womelsdorf. Please visit https://www.biorxiv.org to read.

The script 1_EXAMPLE_SCRIPT_OF_MAGE_USE.m calls several of the main functions used by MAGE, including mage_decomp1.m and mage_filter1.m. These are likely the main functions you are interested in. 

Other functions you may find useful are mage_params2signal.m (converts Nx5 Gabor atom parameter array into a 1-D complex-valued signal that is the weighted sum of the Gabor atoms specified. mage_params2rsignal.m outputs a real-valued synthetic signal.

Three other functions that may prove useful:
ip2g.m outputs the inner product between 2 Gabor atoms, while const_aip_shell_probe_gabors.m outputs a set of Gabor atom parameters that all have a inner product magnitude of aip, given a focus Gabor atom. To visualize these atoms, use make_chirplet.m.
