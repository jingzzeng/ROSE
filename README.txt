=============================================
============== Introduction =================
=============================================

We prepare three demo files, demo1.html, demo2.html, and demo3.html, which illustrate how we implement our proposal and competing methods in an integrated workflow. In these three demo files, we:

  - demo1: In this file, we demonstrate how to implement our proposed ROSE method and reproduce the numerical results in Section 6.1. Additionally, we also demonstrate the codes of other linear dimension reduction competitors, including SEAS-SIR, Lasso-SIR, and Lasso.

  - demo2: In this file, we reproduce the results from three non-linear competing methods in Section 6.2, including autoencoder (Hinton and Salakhutdinov, 2006), deep neural network-based SDR (DNSDR; Chen et al., 2024), and LassoNet (Lemhadri et al., 2021).

  - demo3: In this file, we reproduce the results in Table 2 in Section 7 of the paper, where we analyzed the CCLE data set. For demonstration purpose, we only generate the training-testing data split once and run codes of all competitors once.


======================================================
============== The summary of files ==================
======================================================

- demo1.Rmd, demo2.Rmd, and demo3.Rmd: the Rmd files generating "demo1.html", "demo2.html", and "demo3.html".

- NSDR.py: provided by the authors of Chen et al. (2024). This file contains the main program of implementing the DNSDR approach.

- environment.yml: the YAML file providing the exact configuration of our working Python environment, including environment name, channels, and dependencies.

- data_dictionary.txt: A detailed data dictionary describing all variables in CCLE data and preprocessing steps.


=======================================================
============== The summary of folders ==================
========================================================

* shared_codes: main codes for implementing ROSE and some competitors.

  - msda_utility.R, utility.R: contain utility functions.

  - oracle.R: implements Oracle-ROSE.

  - robustSIR.R: implements ROSE.

  - seas.R: implements SEAS-SIR.


* MLR_codes: main codes for implementing the classic MLR method, two high-dimensional MLR methods, i.e., FGEM and LPEM, and other utility functions. The codes are kindly provided by the authors of Wang et al. (2024).

  - Reference: Wang, N., Deng, K., Mai, Q., and Zhang, X. (2024). Leveraging independence in high-dimensional mixed linear regression. Biometrics, 80(3):ujae103.


* CCLE_data: contains the data sets used in CCLE data analysis.

  - The original files "CCLE_expression.csv" and "CCLE_Drugdata.csv" are provided by the authors of Wang et al. (2024). The file "CCLE_expression.csv" contains the gene expression data of 18,926 genes on 1037 cell lines, and the file "CCLE_Drugdata.csv" contains the response data, describing the information for 504 cell lines.

  - Preprocessing: the original data is preprocessed and the subset with the dimension p = 1000 is saved as "17-AAG_X.rda" and "17-AAG_Y.rda".

  - Data splitting: randomly split the preprocessed data and save the training and testing sets as "dat1.rda".

  - Since DNSDR is implemented in Python environment. For a fair comparison with other competitors implemented in R, we save the reduced training and testing sets produced by the DNSDR method as "dnsdr_17-AAG_train_tf4.csv" and "dnsdr_17-AAG_test_tf4.csv". Then, we load these data into the R environment and make the prediction in R.


* NN_data: contains the training and testing sets generated from demo2.Rmd, which are used for evaluating the performance of three non-linear competing methods, including autoencoder (Hinton and Salakhutdinov, 2006), deep neural network-based SDR (DNSDR; Chen et al., 2024), and LassoNet (Lemhadri et al., 2021).


* demo1_cache, demo2_cache, demo3_cache: cache files produced from the compilation of the Rmd files.

* __pycache__: cache files produced by the NSDR.py file when compiling demo2.Rmd and demo3.Rmd.
