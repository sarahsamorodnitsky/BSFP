% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bsfp.R
\name{bsfp}
\alias{bsfp}
\title{Bayesian Simultaneous Factorization and Prediction (BSFP).}
\usage{
bsfp(
  data,
  Y,
  nninit = TRUE,
  model_params = NULL,
  prior_beta_data_driven = FALSE,
  ranks = NULL,
  scores = NULL,
  nsample,
  progress = TRUE,
  starting_values = NULL,
  save_structures = TRUE,
  save_loadings_scores = TRUE,
  save_predictive_model = TRUE,
  save_imputations = TRUE,
  thinning = 1,
  previous_init = NULL,
  save_init = TRUE,
  save_last_sample = FALSE
)
}
\arguments{
\item{data}{A matrix of lists or a list of matrices that share the same number of
columns. The matrices must be oriented in \eqn{p \times n} orientation. May contain NAs if
there are missing values in the dataset. Datasets not required to have the same number
of features.}

\item{Y}{A matrix of lists or a \eqn{n\times 1} matrix of continuous or binary outcome.
May be \code{NULL} if no outcome is given. May contain NAs if there are missing outcomes.}

\item{nninit}{Boolean determining if nuclear-norm penalized objective is used
to initialize the model. If \code{TRUE}, \code{ranks = NULL}.}

\item{model_params}{A list of hyperparameters for the model.
May be left \code{NULL} if theoretical defaults are desired. Otherwise, must be in
the following named-list form: \code{model_params = list(error_vars = c(), joint_var = double(),
indiv_vars = c(), beta_vars = c(), response_vars = c())}. If we have \eqn{q} sources of data,
\code{error_vars}, \code{indiv_vars} are vectors of length \eqn{q} for each source.
\code{response_vars} must define the shape and rate for the Inverse-Gamma prior of the response variance.
\code{beta_vars} must be of length \eqn{q+1}
to specify the prior variance on the intercept, the prior variance on each joint factor's
contribution to the outcome, and for each individual factor's contribution from each source.
If not specified, must specify prior_beta_data_driven to be TRUE.}

\item{prior_beta_data_driven}{Boolean indicating whether the prior variances on the regression
coefficients should scale with the data. If TRUE, the variance of each regression coefficient is
the ratio of the variance of the outcome over the average variance of each joint or individual factor
(depending on the coefficient). The results from initializing with BIDIFAC will be used, or the user-provided
starting_values will be used. If FALSE, a default of 1 is used. For the intercept, a vague prior
variance of 1e6 is used.}

\item{ranks}{A list of length \eqn{q+1} for the ranks of the joint and individual structures.
Leave \code{NULL} if \code{nninit=TRUE}.}

\item{scores}{Matrix with scores estimated by existing factorization method. Use only
if desired to run Bayesian linear model with scores estimated separately. Otherwise,
leave \code{NULL}.}

\item{nsample}{Integer specifying the number of posterior samples to generate.}

\item{progress}{Boolean indicating if a progress bar be displayed to visualize the progress of the sampler}

\item{starting_values}{List of initial values for Gibbs sampler. If \code{NULL} and \code{nninit=TRUE},
fixes at posterior mode. If \code{NULL} and \code{nninit=FALSE}, simulates from prior distributions.}

\item{save_structures}{Boolean indicating whether the estimated structures should be saved.
Default is TRUE. If TRUE, the estimated overall (joint, individual) structures are saved. Setting to FALSE
may save on memory requirements if inference on the overall structures is not of interest. Note that
the overall structures may be calculated using the estimated loadings and scores, so one may choose to
save those only.}

\item{save_loadings_scores}{Boolean indicating whether estimated loadings and scores should be saved.
Default is TRUE. If TRUE, the estimated (joint, individual) loadings and scores are saved. Setting to FALSE
may save on memory requirements if inference on the loadings and scores is not of interest.}

\item{save_predictive_model}{Boolean indicating whether the estimated regression coefficients, outcome variance
(if outcome is continuous) and underlying continuous variable (if outcome is binary) should be saved. Setting to
FALSE may save on memory requirements if inference on the predictive model is not of interest.}

\item{save_imputations}{Boolean indicating whether the imputed values in the data and/or outcome should be saved.
Setting to FALSE may save on memory requirements if inference on the imputed values is not of interest.}

\item{thinning}{Integer indicating how posterior samples should be thinned. If thinning=10, for example,
every 10th iteration is saved. Default is 1 which implies no thinning. Increasing the thinning may save on memory requirements.}

\item{previous_init}{Path to previous UNIFAC initialization to save on time.}

\item{save_init}{Boolean indicating whether the UNIFAC initialization should be saved. If TRUE, will be saved in current working directory.}

\item{save_last_sample}{Boolean indicating whether the last posterior sample should be saved. This may be useful if more
posterior samples may be needed and this value can be used as in the \code{starting_values} argument.
If \code{save_loadings_scores = TRUE}, this is not necessary.}
}
\value{
Returns a list with the following elements:
\item{data}{Data scaled to error variance 1}
\item{Y}{Response vector, if provided}
\item{J.draw}{List of posterior samples for the estimated joint structure for each source}
\item{A.draw}{List of posterior samples for the estimated individual structure for each source}
\item{S.draw}{List of posterior samples for the overall (joint + individual) structure for each source}
\item{EY.draw}{List of posterior samples for the E(Y|X), i.e. \eqn{\beta_0 + \mathbf{V}\boldsymbol{\beta}_{joint} + \sum_{s=1}^q \mathbf{V}_s \boldsymbol{\beta}_s} for each Gibbs sampling iteration.}
\item{V.draw}{List of posterior samples for joint scores, \eqn{\mathbf{V}}}
\item{U.draw}{List of posterior samples for joint loadings for each source, \eqn{\mathbf{U}_s} for \eqn{s=1,\dots,q}}
\item{W.draw}{List of posterior samples for individual loadings for each source,  \eqn{\mathbf{W}_s} for \eqn{s=1,\dots,q}}
\item{Vs.draw}{List of posterior samples for individual scores for each source, \eqn{\mathbf{V}_s} for \eqn{s=1,\dots,q}}
\item{Xm.draw}{List of predicted values for missing observations in each source \eqn{\mathbf{X}_s} for \eqn{s=1,\dots,q}}
\item{Ym.draw}{List of predicted values for missing outcomes}
\item{Z.draw}{List of draws for latent continuous variable to facilitate Gibbs sampling if outcome is binary}
\item{scores}{Estimated scores provided by a different factorization method in order to run the predictive model}
\item{ranks}{Vector with the estimated joint and individual ranks. \code{ranks[1]} is the estimated joint rank. \code{ranks[2:(q+1)]} correspond to the individual ranks for each source.}
\item{model_params}{List of hyperparameters used in model fitting. If not specified by user, these are the theoretical defaults. If specified by user, returns what was given.}
\item{tau2.draw}{List of posterior samples for the response variance if the response was continuous}
\item{beta.draw}{List of posterior samples for the regression coefficients used in the predictive model}
\item{last.iter}{Last posterior sample for each estimated parameter to use as a future starting value if needed.}
}
\description{
Given multiple sources of data and a continuous or binary outcome measured on \eqn{n} samples,
BSFP can decompose variation across the sources into joint and individual structures.
BSFP simultaneously uses the estimated factors driving these structures to predict an outcome.
BSFP estimates the full posterior distributions of the factors and the predictive model.
Use this function to simulate samples from posterior distributions
of joint and individual structures. This function
can be used in several ways:
(1) Initialize at the theoretical mode of the decomposition. Priors are fixed
at theoretical values. Ranks are determined at initialization. Can include an
outcome or not. If no outcome is included, no burn-in is used. If an outcome
is included, a burn-in must be specified or the default is \code{nsample/2}.
(2) User specifies ranks to estimate. Model is initialized using priors and
then function samples from posterior distributions. User must specify the
hyperparameters for the prior distributions on the structures.
}
\details{
BSFP assumes the features (rows) of each source are centered. It does not require features
to be scaled to overall standard deviation 1. If initializing using the nuclear-norm penalized objective,
the sources will be scaled to have approximately an error variance of 1. Scaling to overall standard
deviation 1 may be overly-conservative in identifying factors.
}
\examples{
# Setting up the data
n <- 50
p.vec <- c(75, 100)
q <- 2

# Setting up the model parameters
true_params <- list(error_vars = c(1,1), # Length must be q
joint_var = 1,
indiv_vars = c(1,1), # Length must be q
beta_vars = c(1, 1, rep(1, q)), # Length must be q+2 (intercept, joint, individual)
response_vars = c(shape = 1, rate = 1))

# Choose ranks
r <- 3
r.vec <- c(3, 3)
ranks <- c(r, r.vec)

# Number of posterior sampling iterations
nsample <- 1000
burnin <- nsample/2
iters_burnin <- (burnin+1):nsample

# Generate data
data.c1 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE)

# Run BSFP for 1000 iterations
bsfp.c1 <- bsfp(data = data.c1$data, Y = data.c1$Y, nsample = nsample)
}
