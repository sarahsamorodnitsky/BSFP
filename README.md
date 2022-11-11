# BSFP - Bayesian Simultaneous Factorization and Prediction

The BSFP R package allows users to integrate multiple sources of omics datasets by decomposing variation into joint and individual structures within a Bayesian framework that allows for full posterior inference. Users may also supply a continuous or binary outcome and use the joint and individual factors driving the estimated structures in a predictive model. The omics sources and outcome vectors may contain missingness, in which case BSFP uses the posterior predictive distribution to perform multiple imputation. The function outputs posterior summaries of the estimated factors and imputed values. 

# Data Preparation

Before applying BSFP to omics data, we recommend centering all features to have mean 0. Users are not required to standardize the features to error variance 1 as doing so will be overly conservative in estimating joint and individual factors. 

Sources should be oriented in a column-linked fashion: the sources should have the same number of columns but not necessarily the same number of rows. 

# Running BSFP

Users are only required to provide the omics sources and an outcome if desired. The BSFP function initialize the model at the posterior mode of the decomposition yielding the estimated number of ranks. The function will return the scaled data, the outcome, estimated ranks, and posterior summaries of the estimated factors after alignment. 

We recommend running BSFP on a high-performance computing system or computing cluster with large memory capacities. 