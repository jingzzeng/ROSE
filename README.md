# ROSE codes

## Introduction

This repository contains the demonstration codes for the paper "Robust Sliced Inverse Regression: Optimal Estimation for Heavy-Tailed Data in High Dimensions". We prepare three demo files, demo1.html, demo2.html, and demo3.html, which illustrate how we implement our proposal and competing methods in an integrated workflow.

  - demo1: Demonstrates how to implement our proposed ROSE method and reproduce the numerical results in Section 6.1. Additionally, we also demonstrate the codes of other linear dimension reduction competitors, including SEAS-SIR, Lasso-SIR, and Lasso.

  - demo2: Reproduces the results from three non-linear competing methods in Section 6.2, including autoencoder (Hinton and Salakhutdinov, 2006), deep neural network-based SDR (DNSDR; Chen et al., 2024), and LassoNet (Lemhadri et al., 2021).

  - demo3: Reproduces the results in Table 2 in Section 7 of the paper, where we analyzed the CCLE data set. For demonstration purpose, we only generate the training-testing data split once and run codes of all competitors once.

## Summary of files

- NSDR.py: provided by the authors of Chen et al. (2024). This file contains the main program of implementing the DNSDR approach.

- environment.yml: the YAML file providing the exact configuration of our working Python environment, including environment name, channels, and dependencies.

- data_dictionary.txt: A detailed data dictionary describing all variables in CCLE data and preprocessing steps.

## Summary of folders

* shared_codes: main codes for implementing ROSE and some competitors.

  - msda_utility.R, utility.R: contain utility functions.

  - oracle.R: implements Oracle-ROSE.

  - robustSIR.R: implements ROSE.

  - seas.R: implements SEAS-SIR.

* MLR_codes: main codes for implementing the classic MLR method, two high-dimensional MLR methods, i.e., FGEM and LPEM, and other utility functions. The codes are kindly provided by the authors of Wang et al. (2024).

* CCLE_data: contains the data sets used in CCLE data analysis. The original files "CCLE_expression.csv" and "CCLE_Drugdata.csv" are provided by the authors of Wang et al. (2024).
  
  - CCLE_expression.csv: the gene expression data of 18,926 genes on 1037 cell lines.
    
  - CCLE_Drugdata.csv: the response data, describing the information for 504 cell lines.

* NN_data: contains the training and testing sets generated from demo2, which are used for evaluating the performance of three non-linear competing methods, including autoencoder (Hinton and Salakhutdinov, 2006), deep neural network-based SDR (DNSDR; Chen et al., 2024), and LassoNet (Lemhadri et al., 2021).

##  Reference

- Zeng, J., Min K., Mai, Q. (2025+). Robust Sliced Inverse Regression: Optimal Estimation for Heavy-Tailed Data in High Dimensions.
  
- Wang, N., Deng, K., Mai, Q., and Zhang, X. (2024). Leveraging independence in high-dimensional mixed linear regression. Biometrics, 80(3):ujae103.

- Hinton, G. E. and Salakhutdinov, R. R. (2006). Reducing the dimensionality of data with neural networks. science, 313(5786):504–507.

- Chen, Y., Jiao, Y., Qiu, R., and Yu, Z. (2024). Deep nonlinear sufficient dimension reduction. The Annals of Statistics, 52(3):1201 – 1226.

- Lemhadri, I., Ruan, F., Abraham, L., and Tibshirani, R. (2021). Lassonet: A neural network with feature sparsity. Journal of Machine Learning Research, 22(127):1–29.
