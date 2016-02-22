This folder contains all the source code for the paper "A Simple Approach to Sparse Clustering" by Ery Arias-Castro and Xiao Pu, and the microarray data sets downloaded from "http://www.stat.cmu.edu/~jiashun/Research/software/GenomicsData/".

Figure1.R: It runs convergence test on Algorithm 1 using syntetic data generated from three mixture of sparse gaussians, and generates Figure1 in the paper.
Figure2.R: It runs simulation to generate Figure2 to illustrate the effectiveness of choosing $s$ using gap statistics.
Figure3.R: It generates Figure3 to give a typical example of the weights that sparse K-means returns.
Figure4.R: It generates Figure 4, which uses a typical dataset (stored in "X.txt") from section 4.3 in the paper and projects the data onto the first two principal components of the submatrix and the whole matrix. This figure illustrates how well the clustering results by SAS_gs and SAS_gss compared with those by sparse Kmeans and IF-PCA.
X.txt: It is a typical dataset generated under the regime in section 4.3, which is used to generate Figure4.
hill_climb.R: It implements our main algorithm, SAS clustering with grid search.
hill_climb_gss: It implements the algorithm SAS clustering with golden section search.
Table1.R: It generates Table 1 to compare hill_climb, hill_climb_GSS, sparse Kmeans and IFPCA under the regime in section 4.1.
Table1_IFPCA.m: This is a Matlab code which performs sparse clustering using IFPCA on the datasets in section 4.1. (We obtained the codes from J. Jin and W. Wang (2014), which can also be downloaded from http://www.stat.cmu.edu/~jiashun/Research/software/HCClustering.)
Table2.R: It generates Table 2 to compare hill_climb, hill_climb_GSS, sparse Kmeans and IFPCA under the regime in section 4.2.
Table2_IFPCA.m: This is a Matlab code which performs sparse clustering using IFPCA on the datasets in section 4.2.
Table3.R: It generates Table 3 to compare hill_climb, hill_climb_GSS, sparse Kmeans and IFPCA under the regime in section 4.3.
Table3_IFPCA.m: This is a Matlab code which performs sparse clustering using IFPCA on the datasets in section 4.3.

Table4_Non_Euclidean: This is a folder which contains the source code for section 4.4. Basically we modified "hill_clim.R" and Witten & Tibshirani's sparseKmeans in the R package "sparcl" , so that we could compare our SAS algorithm with sparseKmeans on high dimensional categorical data.
Table4_Non_Euclidean/hill_climb_pam.R: This file implements our SAS algorithm on categorical data.
Table4_Non_Euclidean/PamSparseCluster.R, Table4_Non_Euclidean/PamSparseClustering.R and Table4_Non_Euclidean/PamSparseClusterPermute: These 3 files modify the sparseKmeans function in the R package "sparcl" , so that we could compare our SAS algorithm with sparseKmeans on high dimensional categorical data.
Table4_Non_Euclidean/Table4: This file generates data under the regime in 4.4, and computes the rand indexes resulted by SAS_gs and sparse-kmedoids.

Table5_Microarray: This is a folder which contains the source code "Table5.R" for section 4.5, and all the 10 microarray datasets used in Table 5, which can also be downloaded from "http://www.stat.cmu.edu/~jiashun/Research/software/GenomicsData/".


