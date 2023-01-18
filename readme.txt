This is the code that generated all the figures in the chapter 'Advances in cytometry gating based on statistical distances and dissimilarities' included in the book 'Statistical Methods at the Forefront of Biomedical Advances'.

We provide both R and Python scripts, and all the necessary datasets. The intended order of execution of the scripts is the following:

1. Lines 1 - 98 in plots_figures_barycenters.R
2. computation_of_W2_barycenter.py
3. Lines 100-163 in plots_figures_barycenters.R
4. computing_distance_matrixes_between_clusters.R
5. trimmed_transport_gates.R

plots_figures_barycenters.R and computation_of_W2_barycenter.py provide computation of Wasserstein barycenters for 2-d cytometry datasets.

computing_distance_matrixes_between_clusters.R contains the computation of different distance matrices between two different cytometries seen as a combination of clusters.

trimmed_transport_gates.R provides optimal transportation of 1d gates.
