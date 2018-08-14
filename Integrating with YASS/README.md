
# Using YASS with GMMSort feature extraction step

The following script adapts the feature extraction step described in [1] to the Yass spike sorting algorithm [2]. 

Briefly, it replaces the PCA step by a wavelet decomposition followed by a weighted-PCA. The weigths are defined using one of the 3 separability metrics discribed in [1].

This script relates to the 'yass.pipeline.run' method. Mainly, it replaces the call for yass.detect.run to yass.detect.wpca.run. To use this pipeline see the following steps:

#### Install Yass-0.9 from https://github.com/paninski-lab/yass. 

#### In the 'yass/detect/' folder add the 'wpca.py' routine from the GMMSort repository: https://github.com/tortlab/GMM-spike-sorting. In the 'yass/threshold' folder add also the 'dimension_reduction_wpca.py'.
    
#### Run the yass calling 'yass.detect.wpca.run' instead of 'yass.detect.run' (see example bellow)

[1] Souza, B.C., Lopes-dos-Santos, V., Bacelo, J. & Tort, A.B.L. (2018). Spike sorting with Gaussian mixture models. bioRxiv.

[2] Lee, J.H., Carlson, D.E., Razaghi, H.S., Yao, W., Goetz, G.A., Hagen, E., ... & Paninski, L. (2017). Yass: Yet another spike sorter. In Advances in Neural Information Processing Systems (pp. 4002-4012).
