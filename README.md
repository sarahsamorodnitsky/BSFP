# BSFP - Bayesian Simultaneous Factorization and Prediction

The BSFP R package allows users to integrate multiple sources of omics datasets by decomposing variation into joint and individual structures within a Bayesian framework that allows for full posterior inference. Users may also supply a continuous or binary outcome and use the joint and individual factors driving the estimated structures in a predictive model. The omics sources and outcome vectors may contain missingness, in which case BSFP uses the posterior predictive distribution to perform multiple imputation. The function outputs posterior summaries of the estimated factors and imputed values.

We recommend running BSFP on a high-performance computing system or computing cluster with large memory capacities. 

To install, run the following:

```
library(devtools)
devtools::install_github("sarahsamorodnitsky/BSFP")
```

This package depends on the following R packages: `MASS`, `Matrix`, `truncnorm`, `infinitefactor`, `svMisc`, `invgamma`, `ggplot2`, `forcats`, and `dplyr`. 

# Examples

Here we will illustrate how to apply the BSFP model to simulated data. The data will be simulated from the assumed model. We will walk through applying BSFP, checking for convergence, aligning the results, and summarizing the results. We will illustrate using BSFP on data with and without missing values. 

Below, we establish some basics of the data we will generate. We will generate $q=2$ sources of data measured on $n=50$ samples. Source 1 will contain $75$ features and source 2 will contain $100$ features. 

We will generate the data to have a true joint rank of $3$ and true individual ranks of $(3,3)$ for each source, respectively. For simplicity, the hyperparameters for the priors of the data-generating model are all fixed at $1$. 

When running BSFP, we will use 1000 posterior sampling iterations with a 500 iteration burn-in. This is likely too low for most real-data applications, but will suffice for this example. 

```{r basic parameters}
# Setting up the data
n <- 50
p.vec <- c(75, 100)
q <- 2

# Setting up the model parameters
true_params <- list(error_vars = c(1,1),
                    joint_var = 1,
                    indiv_vars = c(1,1),
                    beta_vars = c(1, 1, rep(1, q)),
                    response_vars = c(shape = 1, rate = 1))

# Choose ranks
r <- 3
r.vec <- c(3, 3)
ranks <- c(r, r.vec)

# Number of posterior sampling iterations
nsample <- 1000
burnin <- nsample/2
iters_burnin <- (burnin+1):nsample
```

# Setting 1: No Missingness

First we start by generating our 2 sources of data with a continuous response outcome. We do not adjust the signal-to-noise in any way, fixing `s2nX=s2nY=NULL`. Note that in the output of the `bsfp_data` function, the data are stored in a matrix of lists. To access source 1, for example, we would run `data.c1$data[[1,1]]` and source 2 `data.c1$data[[2,1]]`. The response vector, `data.c1$Y`, is stored in a similar way and would be accessed using `data.c1$Y[[1,1]]`. 

```{r setting 1 data generation}
# Generate data
data.c1 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE)

# The simulated observed data is stored in objects data and Y
data.c1$data
data.c1$Y
```

We then run BSFP on the simulated data with the simulated outcome. Note that in general, we do not recommend running BSFP on a local machine. Rather, we recommend using a high-performance computing system as BSFP returns many posterior samples that may require several gigabytes of memory. After running BSFP, we recommend saving the posterior samples and then checking for convergence. 

The function to run BSFP, `bsfp`, has several arguments. In general, users are only required to specify the data (which may be specified in a matrix of lists or using a simple list) and the number of samples. An outcome is not necessarily required if no prediction is desired. The function will, be default, use random matrix theory-motivated hyperparameters in prior distributions. Users may specify their own hyperparameters using the `model_params` argument if they wish. Ranks will be determined by solving a nuclear-norm penalized objective using the `BIDIFAC` function, written by Jun Young Park (2020) [[1]](#1). Alternatively, users may provide their own set of ranks using the `ranks` argument. 

```{r setting 1 run BSFP}
bsfp.c1 <- bsfp(data = data.c1$data, Y = data.c1$Y, nsample = nsample)
```

Once BSFP has finished running, we can check the log-joint density to monitor convergence. If the model has converged, we should see the log-joint density increase and then stabilize. When the model has converged, we should also see that, after a burn-in, the log-joint density values at each posterior sampling iteration appear random and uncorrelated. 

```{r setting 1 convergence}
# Check convergence
log_joint_density_by_iter.c1 <- sapply(1:nsample, function(iter) {
  log_joint_density(data = data.c1$data, Y = data.c1$Y,
                    U.iter = bsfp.c1$U.draw[[iter]],
                    V.iter = bsfp.c1$V.draw[[iter]],
                    W.iter = bsfp.c1$W.draw[[iter]],
                    Vs.iter = bsfp.c1$Vs.draw[[iter]],
                    model_params = bsfp.c1$model_params,
                    ranks = bsfp.c1$ranks,
                    beta.iter = bsfp.c1$beta.draw[[iter]],
                    tau2.iter = bsfp.c1$tau2.draw[[iter]])
})

# Expect to see increase and stabilize over convergence
plot(log_joint_density_by_iter.c1)
```

We can also check for convergence by looking at trace plots of parameters that are NOT affected by rotational, permutation, or sign invariance, namely the joint and individual structures and the conditional expectation of the response vector. Here, we randomly select the (1,1) entry in the joint structure (corresponding to feature 1 and sample 1 from source 1) and plot the posterior samples across Gibbs sampling iterations. For the individual structure, we randomly select source 2 and the (5,10) entry. We should see a "fuzzy caterpillar"-like shape of the trace plots after a burn-in period. 

```{r setting 1 trace plots}
# Check some trace plots of structures

# Source 1 joint structure, (1,1) entry
plot(sapply(1:nsample, function(iter) {
  bsfp.c1$J.draw[[iter]][[1,1]][1,1]
}))

# Source 2 individual structure, (5,10) entry
plot(sapply(1:nsample, function(iter) {
  bsfp.c1$A.draw[[iter]][[2,1]][5,10]
}))
```

We can examine the proportion of variance explained by the estimated structures in the observed datasets using the `var_explained` function.

```{r}
bsfp.c1.var.explained <- var_explained(BSFP.fit = bsfp.c1, iters_burnin = iters_burnin, source.names = c("Expression", "Methylation"))
```

As mentioned above, some model parameters are not identifiable immediately from model fitting due to rotational, permutation, and sign ambiguity. To address this, we modify the MatchAlign algorithm [[2]](#2) to address non-identifiability among the posterior samples. Our code is modified from the R package `infinitefactor` provided by Evan Poworoznek [[3]](#3).

```{r setting 1 alignment}
# Run the alignment algorithm
alignment.c1 <- match_align_bsfp(BSFP.fit = bsfp.c1, y = data.c1$Y,
                                 model_params = bsfp.c1$model_params,
                                 p.vec = p.vec, iters_burnin = iters_burnin)

```

After aligning the factors, we can then summarize them using standard posterior summaries (posterior means, 95\% credible intervals, etc.) We provide the `summarize_factors` function to provide these summaries.


```{r setting 1 posterior summaries}
# Summarize aligned factors
summary.aligned.c1 <- summarize_factors(data = data.c1$data, Y = data.c1$Y,
                                        iters_burnin = iters_burnin,
                                        aligned_results = alignment.c1,
                                        ranks = bsfp.c1$ranks, tau2.draw = bsfp.c1$tau2.draw)
```

We may then want to visualize the posterior summaries of the estimated factors using the `plot_summaries` function. This function requires we specify which structure type (joint or individual) and estimated parameter we'd like to visualize (scores, loadings, regression coefficients). 

```{r setting 1 plotting posterior summaries}
# Plot summaries
plots.joint.scores <- plot_summaries(summary.aligned.c1, structure = "joint", output = "scores")
plots.joint.loadings.source1 <- plot_summaries(summary.aligned.c1, structure = "joint", output = "loadings", source = 1, xlab.name = "Expression")
plots.joint.betas <- plot_summaries(summary.aligned.c1, structure = "joint", output = "betas")
 
plots.individual.scores.source2 <- plot_summaries(summary.aligned.c1, structure = "individual", output = "scores", source = 2)
plots.individual.loadings.source2 <- plot_summaries(summary.aligned.c1, structure = "individual", output = "loadings", source = 2, xlab.name = "Expression")
plots.individuaul.betas.source2 <- plot_summaries(summary.aligned.c1, structure = "individual", output = "betas", source = 2)
 
# View one at a time
plots.joint.scores[[1]] # Joint factor 1
plots.joint.loadings.source1[[2]] # Joint factor 2
plots.joint.betas[[1]] # All regression coefficients for joint factors
 
plots.individual.scores.source2[[1]] # Scores for individual factor 1 from source 2
plots.individual.loadings.source2[[2]] # Loadings for individual factor 2 from source 2
plots.individuaul.betas.source2[[1]] # All regression coefficients for individual factors for source 2
 
# OR, output to a pdf (will output to current working directory)
pdf("Joint_Scores_BSFP.pdf")
plots.joint.scores
dev.off()
```

Finally, we can examine the predictive accuracy of the model fit to the training data. This can also help us gauge convergence of the algorithm. 

```{r setting 1 prediction accuracy}
# Save the fitted values on the training data
Y.pred <- bsfp.c1$EY.draw
Y.pred <- do.call(cbind, lapply(Y.pred, function(i) i[[1,1]]))
 
# Examine convergence
plot(Y.pred[1,])
plot(Y.pred[1,burnin:nsample]) # Add burn-in
 
# Assess prediction accuracy
Y.pred.mean <- rowMeans(Y.pred[,burnin:nsample])
cor(data.c1$Y[[1,1]], Y.pred.mean)
```

# Setting 2: Multiple Imputation

In setting 2, we consider inducing entrywise (missing-at-random) missingness among the data sources to illustrate multiple imputation with BSFP. We start by simulating some data from the model using the same parameters defined previously. We set 10\% of samples in each source to be missing using the `prop_missing` argument. 

```{r setting 2 data generation}
# Generate data
data.c2 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE,
                     missingness = "missingness_in_data", missing_data_type = "entrywise", prop_missing = 0.1)

```

We then fit the BSFP model with the simulated data.

```{r setting 2 run BSFP}
# Run BSFP for 1000 iterations
bsfp.c2 <- bsfp(data = data.c2$missing_data, Y = data.c2$Y, nsample = nsample)
```

We can check convergence using the log-joint density, which will now consider the imputation of missing values in the data.

```{r setting 2 convergence}
# Check convergence
log_joint_density_by_iter.c2 <- sapply(1:nsample, function(iter) {
  log_joint_density(data = data.c2$data, Y = data.c2$Y,
                    U.iter = bsfp.c2$U.draw[[iter]],
                    V.iter = bsfp.c2$V.draw[[iter]],
                    W.iter = bsfp.c2$W.draw[[iter]],
                    Vs.iter = bsfp.c2$Vs.draw[[iter]],
                    model_params = bsfp.c2$model_params,
                    ranks = bsfp.c2$ranks,
                    beta.iter = bsfp.c2$beta.draw[[iter]],
                    tau2.iter = bsfp.c2$tau2.draw[[iter]],
                    Xm.iter = bsfp.c2$Xm.draw[[iter]])
})

# Expect to see increase and stabilize over convergence
plot(log_joint_density_by_iter.c2)
```

We can also check trace plots. We may study the trace plots for the imputed values, too, which are not affected by identifiability issues.

```{r setting 2 trace plots}
# Check some trace plots of structures

# Source 1 joint structure, (1,1) entry
plot(sapply(1:nsample, function(iter) {
  bsfp.c2$J.draw[[iter]][[1,1]][1,1]
}))

# Source 2 individual structure, (5,10) entry
plot(sapply(1:nsample, function(iter) {
  bsfp.c2$A.draw[[iter]][[2,1]][5,10]
}))
```

We can run the alignment algorithm.

```{r setting 2 alignment}
# Run the alignment algorithm
alignment.c2 <- match_align_bsfp(BSFP.fit = bsfp.c2, y = data.c2$Y,
                                 model_params = bsfp.c2$model_params,
                                 p.vec = p.vec, iters_burnin = iters_burnin)
```

Finally, we can study posterior summaries of the estimated model parameters. `summarize_factors` will also output summaries for the posterior samples for imputed values. 

```{r setting 2 posterior summaries}
# Summarize aligned factors
summary.aligned.c2 <- summarize_factors(data = data.c2$missing_data, Y = data.c2$Y,
                                        iters_burnin = iters_burnin,
                                        aligned_results = alignment.c2,
                                        ranks = bsfp.c2$ranks, tau2.draw = bsfp.c2$tau2.draw,
                                        Xm.draw = bsfp.c2$Xm.draw)

```

From here, we may examine the results in a similar manner as Example 1 by looking at posterior summaries for the estimated scores, loadings, and regression coefficients. 

# Example 3: Predicting on Test Data

In this example, we consider fitting BSFP on training data and using the training fit to predict on a held-out test dataset. This involves estimating a new set of joint and individual scores and predicting a previously-unseen response vector. 

We start by generating the data:

```{r}
# Generate data
data.c3 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE)
```

We then split the two sources and the response vector into a training and test dataset.

```{r}
# Split into training and test set
train.c3 <- data.c3$data
train.c3[[1,1]] <- train.c3[[1,1]][,1:(n/2)]
train.c3[[2,1]] <- train.c3[[2,1]][,1:(n/2)]

Y.train.c3 <- data.c3$Y
Y.train.c3[[1,1]] <- Y.train.c3[[1,1]][1:(n/2),,drop=FALSE]


test.c3 <- data.c3$data
test.c3[[1,1]] <- test.c3[[1,1]][,((n/2)+1):n]
test.c3[[2,1]] <- test.c3[[2,1]][,((n/2)+1):n]

Y.test.c3 <- data.c3$Y
Y.test.c3[[1,1]] <- Y.test.c3[[1,1]][((n/2)+1):n,,drop=FALSE]
```

We fit the BSFP model on the training data.

```{r}
# Run BSFP for 1000 iterations
bsfp.train.c3 <- bsfp(data = train.c3, Y = Y.train.c3, nsample = nsample)
```

And predict on the held-out test dataset.

```{r}
# Run BSFP.predict for 1000 iterations on held-out test data
bsfp.test.c3 <- bsfp.predict(bsfp.fit = bsfp.train.c3, test_data = test.c3, Y_test = Y.test.c3,
                             nsample = nsample)
```

We can then check the prediction accuracy.

```{r}
Y.pred <- bsfp.test.c3$EY.draw
Y.pred <- do.call(cbind, lapply(Y.pred, function(i) i[[1,1]]))

# Examine convergence
plot(Y.pred[1,])
plot(Y.pred[1,burnin:nsample]) # Add burn-in

# Assess prediction accuracy
Y.pred.mean <- rowMeans(Y.pred[,burnin:nsample])
cor(Y.test.c3[[1,1]], Y.pred.mean)
```




# References
<a id="1">[1]</a> 
Park, Jun Young, and Eric F. Lock. (2020). 
Integrative factorization of bidimensionally linked matrices. 
Biometrics 76.1 (2020): 61-74.

<a id="2">[2]</a> 
Poworoznek, Evan, Federico Ferrari, and David Dunson. (2021). 
Efficiently resolving rotational ambiguity in Bayesian matrix sampling with matching.
arXiv preprint arXiv:2107.13783.

<a id="3">[3]</a> 
Poworoznek, Evan. (2020). 
Package ‘infinitefactor’. 




