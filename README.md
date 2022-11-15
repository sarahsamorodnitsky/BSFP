# BSFP - Bayesian Simultaneous Factorization and Prediction

The BSFP R package allows users to integrate multiple sources of omics datasets by decomposing variation into joint and individual structures within a Bayesian framework that allows for full posterior inference. Users may also supply a continuous or binary outcome and use the joint and individual factors driving the estimated structures in a predictive model. The omics sources and outcome vectors may contain missingness, in which case BSFP uses the posterior predictive distribution to perform multiple imputation. The function outputs posterior summaries of the estimated factors and imputed values.

We recommend running BSFP on a high-performance computing system or computing cluster with large memory capacities. 

To install, run the following:

```
library(devtools)
devtools::install_github("sarahsamorodnitsky/BSFP")
```





