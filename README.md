# GMM-spike-sorting
The main function 'GMMsort.m' opens a GUI to perform spike sorting and to manually adjust the clusters. 

Alternatively, the function can sort waveforms without opening the GUI (type 'help GMMsort' for help and examples). 
In this case, it uses three other important functions: 

    extract_features.m: performs the feature extraction with wavelet decomposition and weighted-PCA described in the manuscript.
    clusterize.m:       estimates the number and center of clusters and further classify each sample using the extracted features.
    plot_model.m:       plots the resultant classification.

These functions can also be used separately (type 'help extract_features', 'help clusterize' or 'help plot_model' for more information).

The file 'GMMinstructions.pdf' gives instructions for using the GUI.


B. C. Souza

Brain Institute, Natal, Brazil

January, 2018.

Reference: Souza et al., Spike sorting with Gaussian mixture models. Under preparation.
