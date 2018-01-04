# GMM-spike-sorting
Spike sorting with Gaussian mixture models

The main function 'GMMsort.m' can open a GUI to perform spike sorting and to manually adjust the clusters. 

Alternatively, the function can sort waveforms without opening the GUI (type 'help GMMsort' for help and examples). In this case, it uses three other important functions: 
    extract_features.m: perform the feature extraction with wavelet decomposition and weighted-PCA described in the manuscript.
    clusterize.m:       estimate the number center of clusters and classify each sample using the extracted features.
    plot_model.m:       plots the resultant classification.

These functions can also be used separately (type 'help extract_features', 'help clusterize' or 'help plot_model' for more information)

The file 'GMMinstructions.pdf' gives instructions for using the GUI.


B. C. Souza January, 2018
Brain Institute, Natal, Brazil
