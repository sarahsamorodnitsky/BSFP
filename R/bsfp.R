# Bayesian Simultaneous Factorization and Prediction Functions

#' Bayesian Simultaneous Factorization and Prediction (BSFP).
#'
#' Given multiple sources of data and a continuous or binary outcome measured on \eqn{n} samples,
#' BSFP can decompose variation across the sources into joint and individual structures.
#' BSFP simultaneously uses the estimated factors driving these structures to predict an outcome.
#' BSFP estimates the full posterior distributions of the factors and the predictive model.
#' Use this function to simulate samples from posterior distributions
#' of joint and individual structures. This function
#' can be used in several ways:
#' (1) Initialize at the theoretical mode of the decomposition. Priors are fixed
#' at theoretical values. Ranks are determined at initialization. Can include an
#' outcome or not. If no outcome is included, no burn-in is used. If an outcome
#' is included, a burn-in must be specified or the default is \code{nsample/2}.
#' (2) User specifies ranks to estimate. Model is initialized using priors and
#' then function samples from posterior distributions. User must specify the
#' hyperparameters for the prior distributions on the structures.
#' @param data A matrix of lists or a list of matrices that share the same number of
#' columns. The matrices must be oriented in \eqn{p \times n} orientation. May contain NAs if
#' there are missing values in the dataset. Datasets not required to have the same number
#' of features.
#' @param Y A matrix of lists or a \eqn{n\times 1} matrix of continuous or binary outcome.
#' May be \code{NULL} if no outcome is given. May contain NAs if there are missing outcomes.
#' @param nninit Boolean determining if nuclear-norm penalized objective is used
#' to initialize the model. If \code{TRUE}, \code{ranks = NULL}.
#' @param model_params A list of hyperparameters for the model.
#' May be left \code{NULL} if theoretical defaults are desired. Otherwise, must be in
#' the following named-list form: \code{model_params = list(error_vars = c(), joint_var = double(),
#' indiv_vars = c(), beta_vars = c(), response_vars = c())}. If we have \eqn{q} sources of data,
#' \code{error_vars}, \code{indiv_vars} are vectors of length \eqn{q} for each source.
#' \code{response_vars} must define the shape and rate for the Inverse-Gamma prior of the response variance.
#' \code{beta_vars} must be of length \eqn{q+1}
#' to specify the prior variance on the intercept, the prior variance on each joint factor's
#' contribution to the outcome, and for each individual factor's contribution from each source.
#' If not specified, must specify prior_beta_data_driven to be TRUE.
#' @param prior_beta_data_driven Boolean indicating whether the prior variances on the regression
#' coefficients should scale with the data. If TRUE, the variance of each regression coefficient is
#' the ratio of the variance of the outcome over the average variance of each joint or individual factor
#' (depending on the coefficient). The results from initializing with BIDIFAC will be used, or the user-provided
#' starting_values will be used. If FALSE, a default of 1 is used. For the intercept, a vague prior
#' variance of 1e6 is used.
#' @param ranks A list of length \eqn{q+1} for the ranks of the joint and individual structures.
#' Leave \code{NULL} if \code{nninit=TRUE}.
#' @param scores Matrix with scores estimated by existing factorization method. Use only
#' if desired to run Bayesian linear model with scores estimated separately. Otherwise,
#' leave \code{NULL}.
#' @param nsample Integer specifying the number of posterior samples to generate.
#' @param progress Boolean indicating if a progress bar be displayed to visualize the progress of the sampler
#' @param starting_values List of initial values for Gibbs sampler. If \code{NULL} and \code{nninit=TRUE},
#' fixes at posterior mode. If \code{NULL} and \code{nninit=FALSE}, simulates from prior distributions.
#' @param save_structures Boolean indicating whether the estimated structures should be saved.
#' Default is TRUE. If TRUE, the estimated overall (joint, individual) structures are saved. Setting to FALSE
#' may save on memory requirements if inference on the overall structures is not of interest. Note that
#' the overall structures may be calculated using the estimated loadings and scores, so one may choose to
#' save those only.
#' @param save_loadings_scores Boolean indicating whether estimated loadings and scores should be saved.
#' Default is TRUE. If TRUE, the estimated (joint, individual) loadings and scores are saved. Setting to FALSE
#' may save on memory requirements if inference on the loadings and scores is not of interest.
#' @param save_predictive_model Boolean indicating whether the estimated regression coefficients, outcome variance
#' (if outcome is continuous) and underlying continuous variable (if outcome is binary) should be saved. Setting to
#' FALSE may save on memory requirements if inference on the predictive model is not of interest.
#' @param save_imputations Boolean indicating whether the imputed values in the data and/or outcome should be saved.
#' Setting to FALSE may save on memory requirements if inference on the imputed values is not of interest.
#' @param thinning Integer indicating how posterior samples should be thinned. If thinning=10, for example,
#' every 10th iteration is saved. Default is 1 which implies no thinning. Increasing the thinning may save on memory requirements.
#' @param previous_init Path to previous UNIFAC initialization to save on time.
#' @param save_init Boolean indicating whether the UNIFAC initialization should be saved. If TRUE, will be saved in current working directory.
#' @param save_last_sample Boolean indicating whether the last posterior sample should be saved. This may be useful if more
#' posterior samples may be needed and this value can be used as in the \code{starting_values} argument.
#' If \code{save_loadings_scores = TRUE}, this is not necessary.
#'
#' @details BSFP assumes the features (rows) of each source are centered. It does not require features
#' to be scaled to overall standard deviation 1. If initializing using the nuclear-norm penalized objective,
#' the sources will be scaled to have approximately an error variance of 1. Scaling to overall standard
#' deviation 1 may be overly-conservative in identifying factors.
#'
#' @return Returns a list with the following elements:
#' \item{data}{Data scaled to error variance 1}
#' \item{Y}{Response vector, if provided}
#' \item{J.draw}{List of posterior samples for the estimated joint structure for each source}
#' \item{A.draw}{List of posterior samples for the estimated individual structure for each source}
#' \item{S.draw}{List of posterior samples for the overall (joint + individual) structure for each source}
#' \item{EY.draw}{List of posterior samples for the E(Y|X), i.e. \eqn{\beta_0 + \mathbf{V}\boldsymbol{\beta}_{joint} + \sum_{s=1}^q \mathbf{V}_s \boldsymbol{\beta}_s} for each Gibbs sampling iteration.}
#' \item{V.draw}{List of posterior samples for joint scores, \eqn{\mathbf{V}}}
#' \item{U.draw}{List of posterior samples for joint loadings for each source, \eqn{\mathbf{U}_s} for \eqn{s=1,\dots,q}}
#' \item{W.draw}{List of posterior samples for individual loadings for each source,  \eqn{\mathbf{W}_s} for \eqn{s=1,\dots,q}}
#' \item{Vs.draw}{List of posterior samples for individual scores for each source, \eqn{\mathbf{V}_s} for \eqn{s=1,\dots,q}}
#' \item{Xm.draw}{List of predicted values for missing observations in each source \eqn{\mathbf{X}_s} for \eqn{s=1,\dots,q}}
#' \item{Ym.draw}{List of predicted values for missing outcomes}
#' \item{Z.draw}{List of draws for latent continuous variable to facilitate Gibbs sampling if outcome is binary}
#' \item{scores}{Estimated scores provided by a different factorization method in order to run the predictive model}
#' \item{ranks}{Vector with the estimated joint and individual ranks. \code{ranks[1]} is the estimated joint rank. \code{ranks[2:(q+1)]} correspond to the individual ranks for each source.}
#' \item{model_params}{List of hyperparameters used in model fitting. If not specified by user, these are the theoretical defaults. If specified by user, returns what was given.}
#' \item{tau2.draw}{List of posterior samples for the response variance if the response was continuous}
#' \item{beta.draw}{List of posterior samples for the regression coefficients used in the predictive model}
#' \item{last.iter}{Last posterior sample for each estimated parameter to use as a future starting value if needed.}
#'
#' @export
#'
#' @examples
#' # Setting up the data
#' n <- 50
#' p.vec <- c(75, 100)
#' q <- 2
#'
#' # Setting up the model parameters
#' true_params <- list(error_vars = c(1,1), # Length must be q
#' joint_var = 1,
#' indiv_vars = c(1,1), # Length must be q
#' beta_vars = c(1, 1, rep(1, q)), # Length must be q+2 (intercept, joint, individual)
#' response_vars = c(shape = 1, rate = 1))
#'
#' # Choose ranks
#' r <- 3
#' r.vec <- c(3, 3)
#' ranks <- c(r, r.vec)
#'
#' # Number of posterior sampling iterations
#' nsample <- 1000
#' burnin <- nsample/2
#' iters_burnin <- (burnin+1):nsample
#'
#' # Generate data
#' data.c1 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE)
#'
#' # Run BSFP for 1000 iterations
#' bsfp.c1 <- bsfp(data = data.c1$data, Y = data.c1$Y, nsample = nsample)

bsfp <- function(data, Y, nninit = TRUE, model_params = NULL, prior_beta_data_driven = FALSE,
                 ranks = NULL, scores = NULL, nsample, progress = TRUE, starting_values = NULL,
                 save_structures = TRUE, save_loadings_scores = TRUE, save_predictive_model = TRUE,
                 save_imputations = TRUE, thinning = 1, previous_init = NULL, save_init = TRUE, save_last_sample = FALSE) {

  # ---------------------------------------------------------------------------
  # Determining type of input data
  # ---------------------------------------------------------------------------

  # Was the data input as a list?
  if (!("matrix" %in% class(data))) {

    # Save the number of sources
    q <- length(data)

    # Initialize new data matrix
    new_data <- matrix(list(), nrow = q, ncol = 1)

    # Add in the sources
    for (s in 1:q) {
      new_data[[s,1]] <- data[[s]]
    }

    # Rename the new version of the data
    data <- new_data
  }

  # If a response was given, was it input as a matrix of lists or as a matrix?
  if (!is.null(Y)) {
    if (class(Y[1,1]) != "list") {

      # Create a new version of Y
      new_Y <- matrix(list(), nrow = 1, ncol = 1)

      # Input the response
      new_Y[[1,1]] <- Y

      # Return to the Y variable
      Y <- new_Y
    }
  }

  # ---------------------------------------------------------------------------
  # Extracting the dimensions
  # ---------------------------------------------------------------------------

  q <- nrow(data) # Number of sources
  p.vec <- apply(data, 1, function(source) nrow(source[[1]])) # Number of features per source
  p <- sum(p.vec) # Total number of features
  n <- ncol(data[[1,1]]) # Number of subjects

  # Adjust total number of samples to account for thinning
  nsample.thin <- nsample/thinning

  # ---------------------------------------------------------------------------
  # Is there a response vector?
  # ---------------------------------------------------------------------------

  response_given <- !is.null(Y[[1,1]])

  # If so, what kind of response is it?
  if (response_given) {
    Y <- matrix(unlist(Y))

    response_type <- if (all(unique(Y) %in% c(0, 1, NA))) "binary" else "continuous"

    # If there is a response, is there missingness in the outcome?
    missingness_in_response <- any(is.na(Y))

    # Which entries are missing?
    missing_obs_Y <- which(is.na(Y))
  }

  # ---------------------------------------------------------------------------
  # Extracting the model parameters
  # ---------------------------------------------------------------------------

  # If no model parameters are given
  if (is.null(model_params)) {

    # Error variances
    error_vars <- sapply(1:q, function(s) 1)

    # Variance of joint structure
    lambda_joint <- sqrt(sum(p.vec)) + sqrt(n)
    sigma2_joint <- joint_var <- 1/(lambda_joint)

    # Variance of individual structure
    lambda_indiv <- sapply(1:q, function(s) sqrt(p.vec[s]) + sqrt(n))
    sigma2_indiv <- indiv_vars <- 1/(lambda_indiv)

    # For the regression coefficients, beta
    if (prior_beta_data_driven & !nninit) stop("If not initializing with BIDIFAC, must specify hyperparameters for betas.")

    if (!prior_beta_data_driven) {
      lambda2_intercept <- 1e6
      lambda2_joint <- 1
      lambda2_indiv1 <- 1
      lambda2_indiv2 <- 1
      beta_vars <- c(lambda2_intercept, lambda2_joint, lambda2_indiv1, lambda2_indiv2)
    }

    # For the response vector
    shape <- 0.01
    rate <- 0.01
  }

  # If model parameters are given
  if (!is.null(model_params)) {
    error_vars <- model_params$error_vars # Error variances
    sigma2_joint <- joint_var <- model_params$joint_var # Variance of joint structure
    sigma2_indiv <- indiv_vars <- model_params$indiv_vars # Variances of individual structure

    # If a prior on the regression betas is given
    if (!is.null(model_params$beta_vars)) {
      beta_vars <- model_params$beta_vars # Variances on betas
    }

    # If a prior on the response variance is given
    if (!is.null(model_params$response_vars)) {
      response_vars <- model_params$response_vars; shape <- response_vars[1]; rate <- response_vars[2] # Hyperparameters of variance of response
    }

    # If there is a response and model_params is not NULL, there must be beta priors and response variance priors
    if ((!is.null(Y[[1,1]]) & is.null(model_params$beta_vars) & !prior_beta_data_driven)) {
      stop("A response vector is given so please provide hyperparameters for the regression coefficients in the form
           of model_params$beta_vars = c(intercept_prior_var, prior_var_on_joint_factors, prior_var_on_individual_factors_source_1, ..., prior_var_on_individual_factors_source_q")
    }

    # If there is a CONTINUOUS response and model_params is not NULL, there must be beta priors and response variance priors
    if ((!is.null(Y[[1,1]]) & is.null(model_params$response_vars))) {
      if (response_type == "continuous") {
        stop("A response vector is given so please provide hyperparameters for the response variance prior in the form
             of a shape and a rate parameter, i.e. model_params$response_vars = c(shape = a, rate = b).")
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Check for missingness in data
  # ---------------------------------------------------------------------------

  # Check for missingness
  missingness_in_data <- any(sapply(data[,1], function(source) any(is.na(source))))

  # Which entries are missing?
  missing_obs <- lapply(data[,1], function(source) which(is.na(source)))

  # ---------------------------------------------------------------------------
  # Obtaining the ranks
  # ---------------------------------------------------------------------------

  if (nninit) {

    if (is.null(previous_init)) {
      # If there is not missing data
      if (!missingness_in_data) {
        rank_init <- BIDIFAC(data, rmt = TRUE, pbar = FALSE, scale_back = FALSE)
      }

      # If there is missing data
      if (missingness_in_data) {
        rank_init <- impute.BIDIFAC(data = data, rmt = TRUE, pbar = FALSE, scale_back = FALSE)
      }

      # Save the initialization
      if (save_init) save(rank_init, file = paste0("BIDIFAC_initialization_", lubridate::today(), ".rda"))
    }

    if (!is.null(previous_init)) {
      load(previous_init)
    }

    # Print when finished
    print("Posterior mode obtained: joint and individual ranks determined.")

    # Saving the results
    sigma.mat <- rank_init$sigma.mat
    C <- rank_init$C
    r <- Matrix::rankMatrix(C[[1,1]])[[1]] # Joint rank
    I <- rank_init$I
    r.vec <- sapply(1:q, function(s) Matrix::rankMatrix(I[[s,1]])) # Individual ranks

    # Scaling the data
    for (s in 1:q) {
      data[[s,1]] <- data[[s,1]]/sigma.mat[s,]
    }

    # Saving the rank vector
    ranks <- c(r, r.vec)
  }

  if (!nninit) {
    r <- ranks[1]
    r.vec <- ranks[-1]
    sigma.mat <- matrix(nrow = q, 1)
    for (s in 1:q) {
      sigma.mat[s,1] <- 1
    }
  }

  r_total <- r + sum(r.vec)
  n_beta <- 1 + r_total

  # Save the indices for the factors from each source
  rank.inds <- lapply(1:(q+1), function(s) {

    # For any structure,
    if (ranks[s] == 0) {
      NULL
    }

    # For joint structure
    else if (s == 1) {
      if (ranks[s] > 0) {
        1:ranks[s]
      }
    }

    # For each individual structure
    else if (s > 1) {
      (sum(ranks[1:(s-1)])+1):sum(ranks[1:s])
    }

  })

  # ---------------------------------------------------------------------------
  # Storing the posterior samples
  # ---------------------------------------------------------------------------

  # Initialize lists to store posterior samples
  V.draw <- U.draw <- Vs.draw <- W.draw <- beta.draw <- tau2.draw <-
    Z.draw <- Ym.draw <- VStar.draw <- Xm.draw <- NULL

  # If storing the structures, initialize full lists for storage
  if (save_loadings_scores) {
    V.draw <- lapply(1:nsample.thin, function(i) matrix(list(), nrow = 1, ncol = 1))
    U.draw <- lapply(1:nsample.thin, function(i) matrix(list(), nrow = q, ncol = 1))
    Vs.draw <- lapply(1:nsample.thin, function(i) matrix(list(), nrow = 1, ncol = q))
    W.draw <- lapply(1:nsample.thin, function(i) matrix(list(), nrow = q, ncol = q))
  }

  # if (!response_given) {
  #   beta.draw <- tau2.draw <- Z.draw <- Ym.draw <- VStar.draw <- lapply(1:nsample.thin, function(i) matrix(list(), nrow = 1, ncol = 1))
  # }

  # if (!missingness_in_data) {
  #   Xm.draw <- lapply(1:nsample.thin, function(i) matrix(list(), nrow = q, ncol = 1))
  # }

  if (response_given) {
    if (save_predictive_model) {
      beta.draw <- Z.draw <- tau2.draw <- VStar.draw <- lapply(1:nsample.thin, function(i) matrix(list(), nrow = 1, ncol = 1))
    }

    if (save_imputations) Ym.draw <- lapply(1:nsample.thin, function(i) matrix(list(), nrow = 1, ncol = 1))
  }

  if (missingness_in_data) {
    if (save_imputations) Xm.draw <- lapply(1:nsample.thin, function(i) matrix(list(), nrow = q, ncol = 1))
  }

  # ---------------------------------------------------------------------------
  # Initialize V, U, V, W and prior on the regression coefficients
  # ---------------------------------------------------------------------------

  # If initializing with nuclear norm, set initialize values to posterior mode
  if (nninit) {

    V0 <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) {
      svd.joint <- svd(rank_init$C[[1,1]])
      V0[[1,1]] <- (svd.joint$v[,1:r, drop = FALSE]) %*% diag(svd.joint$d[1:r], nrow = r)
    }
    if (r == 0) {
      V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
    }

    U0 <- matrix(list(), nrow = q, ncol = 1)
    Vs0 <- matrix(list(), nrow = 1, ncol = q)
    W0 <- matrix(list(), nrow = q, ncol = q)

    for (s in 1:q) {

      # Initialize U
      if (r > 0) {
        U0[[s,1]] <- svd(rank_init$C[[s,1]])$u[,1:r, drop = FALSE]
      }
      if (r == 0) {
        U0[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
      }

      # Initialize W and V
      if (r.vec[s] > 0) {
        svd.indiv.s <- svd(rank_init$I[[s,1]])
        Vs0[[1,s]] <- (svd.indiv.s$v[,1:r.vec[s], drop = FALSE]) %*% diag(svd.indiv.s$d[1:r.vec[s]], nrow = r.vec[s])
        W0[[s,s]] <- svd.indiv.s$u[,1:r.vec[s], drop = FALSE]

        for (ss in 1:q) {
          if (ss != s) {
            if (r.vec[ss] > 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
            }

            if (r.vec[ss] == 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
            }
          }
        }
      }
      if (r.vec[s] == 0) {
        Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
        W0[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)

        for (ss in 1:q) {
          if (ss != s) {
            if (r.vec[ss] > 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
            }

            if (r.vec[ss] == 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
            }
          }
        }
      }

    }
  }

  # If ranks provided, initialize with prior or use given starting values
  if (!nninit) {

    # If no starting values were provided, initialize from priors
    if (is.null(starting_values)) {
      V0 <- matrix(list(), nrow = 1, ncol = 1)
      if (r > 0) {
        V0[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
      }
      if (r == 0) {
        V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
      }

      U0 <- matrix(list(), nrow = q, ncol = 1)
      Vs0 <- matrix(list(), nrow = 1, ncol = q)
      W0 <- matrix(list(), nrow = q, ncol = q)

      for (s in 1:q) {

        # Initialize U
        if (r > 0) {
          U0[[s,1]] <- matrix(rnorm(p.vec[s]*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = p.vec[s], ncol = r)
        }
        if (r == 0) {
          U0[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }

        # Initialize W and V
        if (r.vec[s] > 0) {
          Vs0[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = n, ncol = r.vec[s])
          W0[[s,s]] <- matrix(rnorm(p.vec[s]*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = p.vec[s], ncol = r.vec[s])

          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }

              if (r.vec[ss] == 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
        if (r.vec[s] == 0) {
          Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
          W0[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)

          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }

              if (r.vec[ss] == 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }

      }
    }

    # If starting values were provided, use them as initial values
    if (!is.null(starting_values)) {
      V0 <- starting_values$V
      U0 <- starting_values$U
      W0 <- starting_values$W
      Vs0 <- starting_values$Vs
    }

  }

  # Establish the option for a prior on the regression coefficients that scales with the data
  if (nninit & prior_beta_data_driven) {

    # Establishing a vague prior on the intercept
    lambda2_intercept <- 1e6

    # Establish prior variance on the joint factors
    lambda2_joint <- var(Y[!is.na(Y)])/mean(apply(V0[[1,1]], 2, var))

    # Establish prior variance on the individual factors
    lambda2_indiv <- var(Y[!is.na(Y)])/sapply(1:q, function(s) mean(apply(Vs0[[1,s]], 2, var)))

    # Combine
    beta_vars <- c(lambda2_intercept, lambda2_joint, lambda2_indiv)

  }

  # Complete the model parameters
  model_params <- list(error_vars = error_vars,
                       joint_var = sigma2_joint,
                       indiv_vars = sigma2_indiv,
                       beta_vars = beta_vars,
                       response_vars = c(shape = shape, rate = rate))

  # If a response is given, set up the variance matrix for the prior of the betas using the ranks
  if (response_given) {
    Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
    beta_vars <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))
    diag(Sigma_beta) <- beta_vars
  }

  if (response_given) {
    # Combining the scores together
    V0.star <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) V0.star[[1,1]] <- V0[[1,1]] else V0.star[[1,1]] <- matrix(nrow = n, ncol = r)

    Vs0.star <- Vs0
    for (s in 1:q) {
      if (r.vec[s] > 0) Vs0.star[[1,s]] <- Vs0[[1,s]] else Vs0.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
    }

    VStar0 <- cbind(1, do.call(cbind, V0.star), do.call(cbind, Vs0.star))

    beta0 <- matrix(MASS::mvrnorm(1, mu = c(rep(0, n_beta)), Sigma = Sigma_beta))
    Z0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = 1))
    tau20 <- matrix(1/rgamma(1, shape = shape, rate = rate))

  }

  # If there is missingness in the data, generate starting values for the missing entries
  if (missingness_in_data) {
    Xm0 <- matrix(list(), ncol = 1, nrow = q)
    for (s in 1:q) {
      if (nninit) {
        Xm0[[s,1]] <- rank_init$X[[s,1]][missing_obs[[s]]]
      }

      if (!nninit) {
        Xm0[[s,1]] <- matrix(0, nrow = length(missing_obs[[s]]), ncol = 1)
      }
    }
  }

  # If there is missingness in Y, generate starting values for the missing entries
  if (response_given) {
    if (missingness_in_response) {
      if (response_type == "continuous") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = sqrt(tau20)))[missing_obs_Y,, drop = FALSE]
      }

      if (response_type == "binary") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rbinom(n, size = 1, prob = pnorm(VStar0 %*% beta0)))[missing_obs_Y,, drop = FALSE]
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Storing the initial values
  # ---------------------------------------------------------------------------

  # Save the initial values in temp variables
  V.iter <- V0
  U.iter <- U0
  Vs.iter <- Vs0
  W.iter <- W0

  if (save_loadings_scores) {
    V.draw[[1]] <- V0
    U.draw[[1]] <- U0
    Vs.draw[[1]] <- Vs0
    W.draw[[1]] <- W0
  }

  if (response_given) {

    # Save initial values in temp variables
    beta.iter <- beta0
    Z.iter <- Z0
    tau2.iter <- tau20
    VStar.iter <- VStar0

    # Store in lists
    if (save_predictive_model) {
      beta.draw[[1]][[1,1]] <- beta0
      Z.draw[[1]][[1,1]] <- Z0
      tau2.draw[[1]][[1,1]] <- tau20
      VStar.draw[[1]][[1,1]] <- VStar0
    }

    if (missingness_in_response) {
      Ym.iter <- Ym0
      if (save_imputations) Ym.draw[[1]][[1,1]] <- Ym0
    }

  }

  if (missingness_in_data) {
    Xm.iter <- Xm0
    if (save_imputations) Xm.draw[[1]] <- Xm0
  }

  # ---------------------------------------------------------------------------
  # Computing the inverses that don't change from iteration to iteration
  # ---------------------------------------------------------------------------

  if (!response_given) {
    # Error variance for X.
    SigmaVInv <- diag(rep(1/error_vars, p.vec))
  }

  if (response_given) {
    if (response_type == "binary") {
      # For V - Combined precisions between data and Z
      SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1))

      # For Vs
      SigmaVsInv <- matrix(list(), nrow = q, ncol = q)

      for (s in 1:q) {
        SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1))
      }
    }

    # For beta - Combined precisions between intercept and all betas
    SigmaBetaInv <- solve(Sigma_beta)
  }

  # ---------------------------------------------------------------------------
  # If structure from another method is given, save the scores as VStar
  # ---------------------------------------------------------------------------

  if (!is.null(scores)) {
    VStar <- cbind(1, scores)
  }

  # ---------------------------------------------------------------------------
  # Start Gibbs sampling!
  # ---------------------------------------------------------------------------

  # Initialize counter for saving posterior samples
  ii <- 1

  for (iter in 1:(nsample-1)) {
    if (progress) svMisc::progress(iter/((nsample-1)/100))

    # Update counter
    if (iter %% thinning == 0) ii <- ii + 1

    # ---------------------------------------------------------------------------
    # Storing the current values of the parameters
    # ---------------------------------------------------------------------------

    # if (save_structures | iter == 1) {
    #   V.iter <- V.draw[[iter]]
    #   U.iter <- U.draw[[iter]]
    #   Vs.iter <- Vs.draw[[iter]]
    #   W.iter <- W.draw[[iter]]
    # }

    if (response_given) {
      # if (save_structures | iter == 1) {
        # # The current values of the betas
        # beta.iter <- beta.draw[[iter]][[1,1]]

        # Creating a matrix of the joint and individual effects
        beta_indiv.iter <- matrix(list(), nrow = q, ncol = 1)

        # Breaking beta down into the intercept,
        beta_intercept.iter <- beta.iter[1,, drop = FALSE]

        # Joint effect
        if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)

        # Individual effects
        if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]

        for (s in 1:q) {
          # If there is no individual effect
          if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)

          # If there is an individual effect
          if (r.vec[s] != 0) {
            if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE]
            if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
          }
        }

        # if (response_type == "binary") {
        #   Z.iter <- Z.draw[[iter]][[1,1]]
        # }
        #
        # if (response_type == "continuous") {
        #   tau2.iter <- tau2.draw[[iter]]
        # }
      # }

      if (missingness_in_response) {
        # Save the current imputations for the missing values
        # Ym.iter <- Ym.draw[[iter]][[1,1]]

        # Creating the completed outcome vector
        Y_complete <- Y

        # Filling in the missing entries for R1 and R2.
        Y_complete[missing_obs_Y,] <- Ym.iter
      }

      if (!missingness_in_response) {
        Y_complete <- Y
      }

    }

    if (missingness_in_data) {
      # Creating the completed matrices.
      X_complete <- data

      # Fill in the completed matrices with the imputed values
      for (s in 1:q) {
        # X_complete[[s,1]][missing_obs[[s]]] <- Xm.draw[[iter]][[s,1]]
        X_complete[[s,1]][missing_obs[[s]]] <- Xm.iter[[s,1]]
      }
    }

    if (!missingness_in_data) {
      X_complete <- data
    }

    # -------------------------------------------------------------------------
    # Computing the inverse that changes with tau2
    # -------------------------------------------------------------------------

    if (response_given) {
      if (response_type == "continuous") {
        # For V - Combined error variances between X1, X2, and Y
        SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1/tau2.iter[[1,1]]))

        # For Vs
        SigmaVsInv <- matrix(list(), nrow = q, ncol = q)

        for (s in 1:q) {
          SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1/tau2.iter[[1,1]]))
        }
      }
    }

    # If estimating the underlying structure
    if (is.null(scores)) {
      # -------------------------------------------------------------------------
      # Posterior sample for V
      # -------------------------------------------------------------------------

      if (r > 0) {
        if (!response_given) {
          # Concatenating Ui's together
          U.iter.combined <- do.call(rbind, U.iter)
          tU_Sigma <- crossprod(U.iter.combined, SigmaVInv)

          # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
          tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)

          # The combined centered Xis with the latent response vector
          X.iter <- do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t))
          Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))

          # V.draw[[iter+1]][[1,1]]
          V.iter[[1,1]] <- t(matrix(sapply(1:n, function(i) {
            bv <-  tU_Sigma %*% X.iter[,i]

            Vi <- MASS::mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
            Vi
          }), nrow = r))
        }

        if (response_given) {
          # Concatenating Ui's together
          U.iter.combined <- rbind(do.call(rbind, U.iter), t(beta_joint.iter))

          # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
          tU_Sigma <- crossprod(U.iter.combined, SigmaVInv)
          tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)

          Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))

          if (response_type == "binary") {
            # The combined centered Xis with the latent response vector
            X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                            t(Z.iter - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
          }

          if (response_type == "continuous") {
            # The combined centered Xis with the latent response vector
            X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                            t(Y_complete - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
          }

          # V.draw[[iter+1]][[1,1]]
          V.iter[[1,1]] <- t(matrix(sapply(1:n, function(i) {
            bv <- tU_Sigma %*% X.iter[,i]

            Vi <- MASS::mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
            Vi
          }), nrow = r))
        }

      }

      if (r == 0) {
        # V.draw[[iter+1]][[1,1]] <- matrix(0, nrow = n, ncol = 1)
        V.iter[[1,1]] <- matrix(0, nrow = n, ncol = 1)
      }

      # Updating the value of V
      # V.iter <- V.draw[[iter+1]]
      if (save_loadings_scores & (iter %% thinning == 0)) V.draw[[ii]][[1,1]] <- V.iter[[1,1]]

      # -------------------------------------------------------------------------
      # Posterior sample for Us
      # -------------------------------------------------------------------------

      if (r > 0) {
        for (s in 1:q) {
          Xs.iter <- X_complete[[s,1]] - W.iter[[s,s]] %*% t(Vs.iter[[1,s]])
          Bu <- solve((1/error_vars[s]) * t(V.iter[[1,1]]) %*% V.iter[[1,1]] + (1/sigma2_joint) * diag(r))
          # U.draw[[iter+1]][[s,1]]
          U.iter[[s,1]] <- t(matrix(sapply(1:p.vec[s], function(j) {
            bu <- (1/error_vars[s]) * t(V.iter[[1,1]]) %*% Xs.iter[j, ]

            U1j <- MASS::mvrnorm(1, mu = Bu %*% bu, Sigma = Bu)
            U1j
          }), nrow = r))
        }
      }

      if (r == 0) {
        for (s in 1:q) {
          # U.draw[[iter+1]][[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
          U.iter[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }
      }

      # U.iter <- U.draw[[iter+1]]
      if (save_loadings_scores & (iter %% thinning == 0)) U.draw[[ii]] <- U.iter

      # -------------------------------------------------------------------------
      # Posterior sample for Vs, s=1,...,q
      # -------------------------------------------------------------------------

      if (!response_given) {
        for (s in 1:q) {
          if (r.vec[s] > 0) {
            Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
            Bvs <- solve((1/error_vars[s]) * t(W.iter[[s,s]]) %*% W.iter[[s,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))

            # Vs.draw[[iter+1]][[1,s]]
            Vs.iter[[1,s]] <- t(matrix(sapply(1:n, function(i) {
              bvs <- (1/error_vars[s]) * t(W.iter[[s,s]]) %*% Xs.iter[, i]

              Vsi <- MASS::mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
              Vsi
            }), nrow = r.vec[s]))
          }

          if (r.vec[s] == 0) {
            # Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
            Vs.iter[[1,s]] <- matrix(0, nrow = n, ncol = 1)
          }
        }
      }

      if (response_given) {
        for (s in 1:q) {
          if (r.vec[s] > 0) {
            # Combined Ws and beta
            W.iter.combined <- rbind(W.iter[[s,s]], t(beta_indiv.iter[[s,1]]))

            tW_Sigma <- crossprod(W.iter.combined, SigmaVsInv[[s,s]])
            tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter.combined)

            Bvs <- solve(tW_Sigma_W + (1/indiv_vars[s]) * diag(r.vec[s]))

            if (response_type == "binary") {
              # Combined centered Xs and Z
              Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                               t(Z.iter - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter -
                                   do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
            }

            if (response_type == "continuous") {
              # Combined centered Xs and Y
              Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                               t(Y_complete - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter -
                                   do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
            }

            # Vs.draw[[iter+1]][[1,s]]
            Vs.iter[[1,s]] <- t(matrix(sapply(1:n, function(i) {
              bvs <- tW_Sigma %*% Xs.iter[, i]

              Vsi <- MASS::mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
              Vsi
            }), nrow = r.vec[s]))
          }

          if (r.vec[s] == 0) {
            # Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
            Vs.iter[[1,s]] <- matrix(0, nrow = n, ncol = 1)
          }
        }
      }

      # Update the current value of V
      # Vs.iter <- Vs.draw[[iter+1]]
      if (save_loadings_scores & (iter %% thinning == 0)) Vs.draw[[ii]] <- Vs.iter

      # Combine current values of V and V.
      V.iter.star.joint <- V.iter
      if (r == 0) {
        V.iter.star.joint[[1,1]] <- matrix(nrow = n, ncol = r)
      }

      Vs.iter.star <- Vs.iter
      for (s in 1:q) {
        if (r.vec[s] == 0) {
          Vs.iter.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
        }
      }

      VStar.iter <- cbind(1, do.call(cbind, V.iter.star.joint), do.call(cbind, Vs.iter.star))

      # Save the current VStar
      if (save_loadings_scores & (iter %% thinning == 0)) VStar.draw[[ii]][[1,1]] <- VStar.iter

      # -------------------------------------------------------------------------
      # Posterior sample for W
      # -------------------------------------------------------------------------

      for (s in 1:q) {
        if (r.vec[s] > 0) {
          Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
          Bws <- solve((1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Vs.iter[[1,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))

          # W.draw[[iter+1]][[s,s]]
          W.iter[[s,s]] <- t(matrix(sapply(1:p.vec[s], function(j) {
            bws <- (1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Xs.iter[j,]

            Wsj <- MASS::mvrnorm(1, mu = Bws %*% bws, Sigma = Bws)
            Wsj
          }), nrow = r.vec[s]))

          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                # W.draw[[iter+1]][[s,ss]]
                W.iter[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }

              if (r.vec[ss] == 0) {
                # W.draw[[iter+1]][[s,ss]]
                W.iter[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }

        if (r.vec[s] == 0) {
          # W.draw[[iter+1]][[s,s]]
          W.iter[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)

          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                # W.draw[[iter+1]][[s,ss]]
                W.iter[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }

              if (r.vec[ss] == 0) {
                # W.draw[[iter+1]][[s,ss]]
                W.iter[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
      }

      # Update the current value of W
      if (save_loadings_scores & (iter %% thinning == 0)) W.draw[[ii]] <- W.iter

    }

    # If structure from another method is provided
    if (!is.null(scores)) {
      VStar.iter <- VStar
    }

    # -------------------------------------------------------------------------
    # Posterior sample for beta
    # -------------------------------------------------------------------------

    if (response_given) {

      if (response_type == "binary") {
        Bbeta <- solve(t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
        bbeta <- t(VStar.iter) %*% Z.iter
      }

      if (response_type == "continuous") {
        Bbeta <- solve((1/tau2.iter[[1,1]]) * t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
        bbeta <- (1/tau2.iter[[1,1]]) * t(VStar.iter) %*% Y_complete
      }

      # beta.draw[[iter+1]][[1,1]]
      beta.iter <- matrix(MASS::mvrnorm(1, mu = Bbeta %*% bbeta, Sigma = Bbeta), ncol = 1)

      # Update the current value of beta
      if (save_predictive_model & (iter %% thinning == 0)) beta.draw[[ii]][[1,1]] <- beta.iter

      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), ncol = 1, nrow = q)

      # Breaking beta down into the intercept
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]

      # Joint effect
      if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)

      # Individual effects
      if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]

      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)

        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE]
          if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }
    }

    # -------------------------------------------------------------------------
    # Posterior sample for tau2
    # -------------------------------------------------------------------------

    if (response_given) {
      if (response_type == "continuous") {
        # tau2.draw[[iter+1]][[1,1]]
        tau2.iter[[1,1]] <- matrix(1/rgamma(1, shape = shape + (n/2), rate = rate + 0.5 * sum((Y_complete - VStar.iter %*% beta.iter)^2)))

        # Update the current value of tau2
        if (save_predictive_model & (iter %% thinning == 0)) tau2.draw[[ii]] <- tau2.iter
      }
    }

    # -------------------------------------------------------------------------
    # Posterior sample for latent continuous response Z
    # -------------------------------------------------------------------------

    if (response_given) {
      if (response_type == "binary") {
        # Z.draw[[iter+1]][[1,1]]
        Z.iter <- matrix(sapply(1:n, function(i) {
          if (Y_complete[i,] == 1) {
            truncnorm::rtruncnorm(1, a = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          } else {
            truncnorm::rtruncnorm(1, b = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          }
        }), ncol = 1)

        if (save_predictive_model & (iter %% thinning == 0)) Z.draw[[ii]] <- Z.iter
      }
    }

    # -------------------------------------------------------------------------
    # Impute missing data
    # -------------------------------------------------------------------------

    if (response_given) {
      if (missingness_in_response) {
        if (response_type == "continuous") {
          # Ym.draw[[iter+1]][[1,1]] <- matrix(rnorm(n, mean = VStar.iter %*% beta.iter, sd = sqrt(tau2.iter[[1,1]])), ncol = 1)[missing_obs_Y,, drop = FALSE]
          Ym.iter <- matrix(rnorm(n, mean = VStar.iter %*% beta.iter, sd = sqrt(tau2.iter[[1,1]])), ncol = 1)[missing_obs_Y,, drop = FALSE]
        }

        if (response_type == "binary") {
          # Ym.draw[[iter+1]][[1,1]] <- matrix(rbinom(n, size = 1, prob = pnorm(VStar.iter %*% beta.iter)), ncol = 1)[missing_obs_Y,, drop = FALSE]
          Ym.iter <- matrix(rbinom(n, size = 1, prob = pnorm(VStar.iter %*% beta.iter)), ncol = 1)[missing_obs_Y,, drop = FALSE]
        }

        if (save_imputations & (iter %% thinning == 0)) Ym.draw[[ii]][[1,1]] <- Ym.iter
      }
    }

    if (missingness_in_data) {
      for (s in 1:q) {
        Es <-  matrix(rnorm(p.vec[s]*n, 0, sqrt(error_vars[s])), nrow = p.vec[s], ncol = n)
        # Xm.draw[[iter+1]][[s,1]] <- matrix((U.iter[[s,1]] %*% t(V.iter[[1,1]]) + W.iter[[s,s]] %*% t(Vs.iter[[1,s]]) + Es)[missing_obs[[s]]])
        Xm.iter[[s,1]] <- matrix((U.iter[[s,1]] %*% t(V.iter[[1,1]]) + W.iter[[s,s]] %*% t(Vs.iter[[1,s]]) + Es)[missing_obs[[s]]])
      }

      if (save_imputations & (iter %% thinning == 0)) Xm.draw[[ii]] <- Xm.iter
    }
  }

  # ---------------------------------------------------------------------------
  # Scaling the imputed data back to the original data scale
  # ---------------------------------------------------------------------------

  # If there is any missingness
  if (missingness_in_data) {
    if (save_imputations) {
      for (iter in 1:nsample.thin) {
        for (s in 1:q) {
          Xm.draw[[iter]][[s,1]] <- Xm.draw[[iter]][[s,1]] * sigma.mat[s,1]
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Calculating the joint and individual structure, scaled to the data
  # ---------------------------------------------------------------------------

  # Initialize lists to store the overall structures
  J.draw <- A.draw <- S.draw <- EY.draw <- NULL

  if (save_structures) {
    if (is.null(scores)) {
      # Storing the structures at each Gibbs sampling iteration
      J.draw <- A.draw <- S.draw <- lapply(1:nsample.thin, function(i) matrix(list(), nrow = q, ncol = 1))

      # Calculating the structure for Y at each Gibbs sampling iteration
      EY.draw <- lapply(1:nsample.thin, function(i) matrix(list(), nrow = 1, ncol = 1))

      for (iter in 1:nsample.thin) {
        for (s in 1:q) {
          # Calculating the joint structure and scaling by sigma.mat
          J.draw[[iter]][[s,1]] <- (U.draw[[iter]][[s,1]] %*% t(V.draw[[iter]][[1,1]])) * sigma.mat[s,1]

          # Calculating the individual structure and scaling by sigma.mat
          A.draw[[iter]][[s,1]] <- (W.draw[[iter]][[s,s]] %*% t(Vs.draw[[iter]][[1,s]])) * sigma.mat[s,1]

          # Calculate the overall structure (joint + individual
          S.draw[[iter]][[s,1]] <- J.draw[[iter]][[s,1]] + A.draw[[iter]][[s,1]]
        }

        # Calculate the structure for Y
        if (response_given) {
          if (response_type == "continuous") {
            EY.draw[[iter]][[1,1]] <- VStar.draw[[iter]][[1,1]] %*% beta.draw[[iter]][[1,1]]
          }

          if (response_type == "binary") {
            EY.draw[[iter]][[1,1]] <- pnorm(VStar.draw[[iter]][[1,1]] %*% beta.draw[[iter]][[1,1]])
          }
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Save the last sample if more posterior samples are needed
  # ---------------------------------------------------------------------------

  last.iter <- list(V.iter = NULL, U.iter = NULL, Vs.iter = NULL, W.iter = NULL,
                    beta.iter = NULL, tau2.iter = NULL, Z.iter = NULL,
                    Ym.iter = NULL, Xm.iter = NULL)

  if (save_last_sample) {
    last.iter$V.iter <- V.iter
    last.iter$U.iter <- U.iter
    last.iter$Vs.iter <- Vs.iter
    last.iter$W.iter <- W.iter

    if (response_given) {
      last.iter$beta.iter <- beta.iter

      if (response_type == "continuous") last.iter$tau2.iter <- tau2.iter
      if (response_type == "binary") last.iter$Z.iter <- Z.iter
      if (missingness_in_response) last.iter$Ym.iter <- Ym.iter
    }

    if (missingness_in_data) last.iter$Xm.iter <- Xm.iter
  }

  # Return
  list(data = data, # Returning the scaled version of the data
       Y = Y, # Return the response vector
       sigma.mat = sigma.mat, # Scaling factors
       J.draw = J.draw, A.draw = A.draw, S.draw = S.draw, EY.draw = EY.draw, # Underlying structure
       V.draw = V.draw, U.draw = U.draw, W.draw = W.draw, Vs.draw = Vs.draw, # Components of the structure
       Xm.draw = Xm.draw, Ym.draw = Ym.draw, Z.draw = Z.draw, # Missing data imputation
       scores = scores, # Scores if provided by another method
       ranks = c(r, r.vec), # Ranks
       model_params = model_params, # Parameters used in model
       tau2.draw = tau2.draw, beta.draw = beta.draw,  # Regression parameters
       last.iter = last.iter) # Last posterior sample for each parameter
}

#' bsfp.predict
#'
#' Gibbs sampling algorithm for sampling new scores for unseen test data.
#' Given BSFP was fit on training data (training sources and outcome), sample
#' the scores and regression coefficient for held-out test data sources and outcome.
#' @param bsfp.fit Results from fitting \code{bsfp} on training data.
#' @param test_data Matrix-list dataset of held-out test data.
#' @param Y_test Column vector with outcome on test samples or \code{NULL}
#' @param model_params May be \code{NULL} if \code{model_params=NULL} in \code{bsfp} fit.
#' Otherwise, specify as \code{(error_vars, joint_vars, indiv_vars, beta_vars, response_vars)}.
#' @param nsample Integer specifing number of Gibbs sampling iterations
#' @param progress Boolean determining if progress of the sampler be displayed
#' @param starting_values List of starting values for \eqn{\mathbf{V}, \mathbf{U}_s, \mathbf{W}_s, \mathbf{V}_s} for \eqn{s=1,\dots, q}.
#' If \code{NULL}, initialize from prior.
#'
#' @details Generate new scores for held-out test data based on a
#' training fit of BSFP. Uses the estimated ranks and joint and individual loadings. Cannot
#' be used if missing values are present in test data.
#'
#' @return Returns a list with the following parameters:
#' \item{test_data}{Test data provided by user}
#' \item{Y_test}{Test response provided by user}
#' \item{J.draw}{List of posterior samples for the estimated joint structure for each source}
#' \item{A.draw}{List of posterior samples for the estimated individual structure for each source}
#' \item{S.draw}{List of posterior samples for the overall (joint + individual) structure for each source}
#' \item{EY.draw}{List of posterior samples for the E(Y|X), i.e. \eqn{\beta_0 + \mathbf{V}\boldsymbol{\beta}_{joint} + \sum_{s=1}^q \mathbf{V}_s \boldsymbol{\beta}_s} for each Gibbs sampling iteration.}
#' \item{V.draw}{List of posterior samples for joint scores, \eqn{\mathbf{V}}}
#' \item{U.train}{List of posterior samples for joint loadings for each source, \eqn{\mathbf{U}_s} for \eqn{s=1,\dots,q} given by the training BSFP fit}
#' \item{W.train}{List of posterior samples for individual loadings for each source,  \eqn{\mathbf{W}_s} for \eqn{s=1,\dots,q} given by the training BSFP fit}
#' \item{Vs.draw}{List of posterior samples for individual scores for each source, \eqn{\mathbf{V}_s} for \eqn{s=1,\dots,q}}
#' \item{ranks}{Vector with the estimated joint and individual ranks. \code{ranks[1]} is the estimated joint rank. \code{ranks[2:(q+1)]} correspond to the individual ranks for each source.}
#' \item{tau2.train}{List of posterior samples for the response variance if the response was continuous given by training BSFP fit}
#' \item{beta.train}{List of posterior samples for the regression coefficients used in the predictive model given by training BSFP fit}
#'
#' @export
#'
#' @examples
#' # Setting up the data
#' n <- 100
#' p.vec <- c(75, 100)
#' q <- 2
#'
#' # Setting up the model parameters
#' true_params <- list(error_vars = c(1,1), # Length must be q
#' joint_var = 1,
#' indiv_vars = c(1,1), # Length must be q
#' beta_vars = c(1, 1, rep(1, q)), # Length must be q+2 (intercept, joint, individual)
#' response_vars = c(shape = 1, rate = 1))
#'
#' # Choose ranks
#' r <- 3
#' r.vec <- c(3, 3)
#' ranks <- c(r, r.vec)
#'
#' # Number of posterior sampling iterations
#' nsample <- 1000
#' burnin <- nsample/2
#' iters_burnin <- (burnin+1):nsample
#'
#' # Generate data
#' data.c3 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE)
#'
#' # Split into training and test set
#' train.c3 <- data.c3$data
#' train.c3[[1,1]] <- train.c3[[1,1]][,1:(n/2)]
#' train.c3[[2,1]] <- train.c3[[2,1]][,1:(n/2)]
#'
#' Y.train.c3 <- data.c3$Y
#' Y.train.c3[[1,1]] <- Y.train.c3[[1,1]][1:(n/2),,drop=FALSE]
#'
#' test.c3 <- data.c3$data
#' test.c3[[1,1]] <- test.c3[[1,1]][,((n/2)+1):n]
#' test.c3[[2,1]] <- test.c3[[2,1]][,((n/2)+1):n]
#'
#' Y.test.c3 <- data.c3$Y
#' Y.test.c3[[1,1]] <- Y.test.c3[[1,1]][((n/2)+1):n,,drop=FALSE]
#'
#' # Run BSFP for 1000 iterations
#' bsfp.train.c3 <- bsfp(data = train.c3, Y = Y.train.c3, nsample = nsample)
#'
#' # Run BSFP.predict for 1000 iterations on held-out test data
#' bsfp.test.c3 <- bsfp.predict(bsfp.fit = bsfp.train.c3, test_data = test.c3, Y_test = Y.test.c3, nsample = nsample)

bsfp.predict <- function(bsfp.fit, test_data, Y_test, model_params = NULL, nsample, progress = TRUE, starting_values = NULL) {

  # ---------------------------------------------------------------------------
  # Extracting the dimensions
  # ---------------------------------------------------------------------------

  q <- nrow(test_data) # Number of sources
  p.vec <- apply(test_data, 1, function(source) nrow(source[[1]])) # Number of features per source
  p <- sum(p.vec) # Total number of features
  n <- ncol(test_data[[1,1]]) # Number of subjects

  # ---------------------------------------------------------------------------
  # Extracting the model parameters
  # ---------------------------------------------------------------------------

  if (is.null(model_params)) {
    model_params <- bsfp.fit$model_params
  }

  error_vars <- model_params$error_vars # Error variances
  sigma2_joint <- joint_var <- model_params$joint_var # Variance of joint structure
  sigma2_indiv <- indiv_vars <- model_params$indiv_vars # Variances of individual structure
  beta_vars <- model_params$beta_vars # Variances on betas
  response_vars <- model_params$response_vars; shape <- response_vars[1]; rate <- response_vars[2] # Hyperparameters of variance of response

  # ---------------------------------------------------------------------------
  # Is there a response vector?
  # ---------------------------------------------------------------------------

  response_given <- !is.null(Y_test[[1,1]])

  # If so, what kind of response is it?
  if (response_given) {
    Y_test <- matrix(unlist(Y_test))

    response_type <- if (all(unique(Y_test) %in% c(0, 1, NA))) "binary" else "continuous"
  }

  # ---------------------------------------------------------------------------
  # Obtaining the ranks
  # ---------------------------------------------------------------------------

  # Saving the ranks from the training data fit
  ranks <- bsfp.fit$ranks

  r <- ranks[1]
  r.vec <- ranks[-1]
  sigma.mat <- matrix(nrow = q, 1)
  for (s in 1:q) {
    sigma.mat[s,1] <- 1
  }

  r_total <- r + sum(r.vec)
  n_beta <- 1 + r_total

  # If a response is given, set up the variance matrix for the prior of the betas using the ranks
  if (response_given) {
    Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
    beta_vars <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))
    diag(Sigma_beta) <- beta_vars
  }

  # ---------------------------------------------------------------------------
  # Save the loadings, betas, and outcome variance estimated on the training data
  # ---------------------------------------------------------------------------

  # Save the joint loadings
  U.train <- bsfp.fit$U.draw

  # Save the individual loadings
  W.train <- bsfp.fit$W.draw

  # Save the betas
  beta.train <- bsfp.fit$beta.draw

  # Save the response variance
  tau2.train <- unlist(bsfp.fit$tau2.draw)

  # ---------------------------------------------------------------------------
  # Storing the posterior samples
  # ---------------------------------------------------------------------------

  V.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  Vs.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = q))

  if (!response_given) {
    Z.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }

  if (response_given) {
    Z.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }

  # ---------------------------------------------------------------------------
  # Initialize V, Vs
  # ---------------------------------------------------------------------------

  # If no starting values were provided, initialize from priors
  if (is.null(starting_values)) {
    V0 <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) {
      V0[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
    }
    if (r == 0) {
      V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
    }

    Vs0 <- matrix(list(), nrow = 1, ncol = q)

    for (s in 1:q) {

      # Initialize W and V
      if (r.vec[s] > 0) {
        Vs0[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = n, ncol = r.vec[s])

      }
      if (r.vec[s] == 0) {
        Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)

      }

    }
  }

  # If starting values were provided, use them as initial values
  if (!is.null(starting_values)) {
    V0 <- starting_values$V
    Vs0 <- starting_values$Vs
  }

  if (response_given) {
    # Combining the scores together
    V0.star <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) V0.star[[1,1]] <- V0[[1,1]] else V0.star[[1,1]] <- matrix(nrow = n, ncol = r)

    Vs0.star <- Vs0
    for (s in 1:q) {
      if (r.vec[s] > 0) Vs0.star[[1,s]] <- Vs0[[1,s]] else Vs0.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
    }

    VStar0 <- cbind(1, do.call(cbind, V0.star), do.call(cbind, Vs0.star))
    Z0 <- matrix(rnorm(n, mean = VStar0 %*% beta.train[[1]][[1,1]], sd = 1))

  }

  # ---------------------------------------------------------------------------
  # Storing the initial values
  # ---------------------------------------------------------------------------

  V.draw[[1]] <- V0
  Vs.draw[[1]] <- Vs0

  if (response_given) {
    Z.draw[[1]][[1,1]] <- Z0
    VStar.draw[[1]][[1,1]] <- VStar0
  }

  # ---------------------------------------------------------------------------
  # Computing the inverses that don't change from iteration to iteration
  # ---------------------------------------------------------------------------

  if (!response_given) {
    # Error variance for X.
    SigmaVInv <- diag(rep(1/error_vars, p.vec))
  }

  if (response_given) {
    if (response_type == "binary") {
      # For V - Combined precisions between data and Z
      SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1))

      # For Vs
      SigmaVsInv <- matrix(list(), nrow = q, ncol = q)

      for (s in 1:q) {
        SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1))
      }
    }

    # For beta - Combined precisions between intercept and all betas
    SigmaBetaInv <- solve(Sigma_beta)
  }

  # ---------------------------------------------------------------------------
  # Start Gibbs sampling!
  # ---------------------------------------------------------------------------

  for (iter in 1:(nsample-1)) {
    if (progress) svMisc::progress(iter/((nsample-1)/100))

    # ---------------------------------------------------------------------------
    # Storing the current values of the parameters
    # ---------------------------------------------------------------------------

    V.iter <- V.draw[[iter]]
    U.iter <- U.train[[iter]]
    Vs.iter <- Vs.draw[[iter]]
    W.iter <- W.train[[iter]]

    if (response_given) {
      # The current values of the betas
      beta.iter <- beta.train[[iter]][[1,1]]

      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), nrow = q, ncol = 1)

      # Breaking beta down into the intercept,
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]

      # Joint effect
      if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)

      # Individual effects
      if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]

      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)

        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE]
          if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }

      if (response_type == "binary") {
        Z.iter <- Z.draw[[iter]][[1,1]]
      }

      if (response_type == "continuous") {
        tau2.iter <- tau2.train[iter]
      }

      Y_complete <- Y_test

    }

    X_complete <- test_data

    # -------------------------------------------------------------------------
    # Computing the inverse that changes with tau2
    # -------------------------------------------------------------------------

    if (response_given) {
      if (response_type == "continuous") {
        # For V - Combined error variances between X1, X2, and Y
        SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1/tau2.iter))

        # For Vs
        SigmaVsInv <- matrix(list(), nrow = q, ncol = q)

        for (s in 1:q) {
          SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1/tau2.iter))
        }
      }
    }

    # -------------------------------------------------------------------------
    # Posterior sample for V
    # -------------------------------------------------------------------------

    if (r > 0) {
      if (!response_given) {
        # Concatenating Ui's together
        U.iter.combined <- do.call(rbind, U.iter)
        tU_Sigma <- crossprod(U.iter.combined, SigmaVInv)

        # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
        tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)

        # The combined centered Xis with the latent response vector
        X.iter <- do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t))
        Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))

        V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
          bv <-  tU_Sigma %*% X.iter[,i]

          Vi <- MASS::mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
          Vi
        }), nrow = r))
      }

      if (response_given) {
        # Concatenating Ui's together
        U.iter.combined <- rbind(do.call(rbind, U.iter), t(beta_joint.iter))

        # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
        tU_Sigma <- crossprod(U.iter.combined, SigmaVInv)
        tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)

        Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))

        if (response_type == "binary") {
          # The combined centered Xis with the latent response vector
          X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                          t(Z.iter - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
        }

        if (response_type == "continuous") {
          # The combined centered Xis with the latent response vector
          X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                          t(Y_complete - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
        }

        V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
          bv <- tU_Sigma %*% X.iter[,i]

          Vi <- MASS::mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
          Vi
        }), nrow = r))
      }

    }

    if (r == 0) {
      V.draw[[iter+1]][[1,1]] <- matrix(0, nrow = n, ncol = 1)
    }

    # Updating the value of V
    V.iter <- V.draw[[iter+1]]

    # -------------------------------------------------------------------------
    # Save the training U, joint loadings
    # -------------------------------------------------------------------------

    U.iter <- U.train[[iter+1]]

    # -------------------------------------------------------------------------
    # Posterior sample for Vs, s=1,...,q
    # -------------------------------------------------------------------------

    if (!response_given) {
      for (s in 1:q) {
        if (r.vec[s] > 0) {
          Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
          Bvs <- solve((1/error_vars[s]) * t(W.iter[[s,s]]) %*% W.iter[[s,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))

          Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
            bvs <- (1/error_vars[s]) * t(W.iter[[s,s]]) %*% Xs.iter[, i]

            Vsi <- MASS::mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
            Vsi
          }), nrow = r.vec[s]))
        }

        if (r.vec[s] == 0) {
          Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
        }
      }
    }

    if (response_given) {
      for (s in 1:q) {
        if (r.vec[s] > 0) {
          # Combined Ws and beta
          W.iter.combined <- rbind(W.iter[[s,s]], t(beta_indiv.iter[[s,1]]))

          tW_Sigma <- crossprod(W.iter.combined, SigmaVsInv[[s,s]])
          tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter.combined)

          Bvs <- solve(tW_Sigma_W + (1/indiv_vars[s]) * diag(r.vec[s]))

          if (response_type == "binary") {
            # Combined centered Xs and Z
            Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                             t(Z.iter - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter -
                                 do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
          }

          if (response_type == "continuous") {
            # Combined centered Xs and Y
            Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                             t(Y_complete - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter -
                                 do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
          }

          Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
            bvs <- tW_Sigma %*% Xs.iter[, i]

            Vsi <- MASS::mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
            Vsi
          }), nrow = r.vec[s]))
        }

        if (r.vec[s] == 0) {
          Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
        }
      }
    }

    # Update the current value of V
    Vs.iter <- Vs.draw[[iter+1]]

    # Combine current values of V and V.
    V.iter.star.joint <- V.iter
    if (r == 0) {
      V.iter.star.joint[[1,1]] <- matrix(nrow = n, ncol = r)
    }

    Vs.iter.star <- Vs.iter
    for (s in 1:q) {
      if (r.vec[s] == 0) {
        Vs.iter.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
      }
    }

    VStar.iter <- cbind(1, do.call(cbind, V.iter.star.joint), do.call(cbind, Vs.iter.star))

    # Save the current VStar
    VStar.draw[[iter+1]][[1,1]] <- VStar.iter

    # -------------------------------------------------------------------------
    # Save training sample for W, individual loadings
    # -------------------------------------------------------------------------

    # Update the current value of W
    W.iter <- W.train[[iter+1]]

    # -------------------------------------------------------------------------
    # Posterior sample for beta
    # -------------------------------------------------------------------------

    if (response_given) {

      # Update the current value of beta
      beta.iter <- beta.train[[iter+1]][[1,1]]

      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), ncol = 1, nrow = q)

      # Breaking beta down into the intercept
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]

      # Joint effect
      if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)

      # Individual effects
      if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]

      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)

        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE]
          if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }
    }

    # -------------------------------------------------------------------------
    # Save the training tau2
    # -------------------------------------------------------------------------

    if (response_given) {
      if (response_type == "continuous") {
        # Update the current value of tau2
        tau2.iter <- tau2.train[iter+1]
      }
    }

    # -------------------------------------------------------------------------
    # Posterior sample for latent continuous response Z
    # -------------------------------------------------------------------------

    if (response_given) {
      if (response_type == "binary") {
        Z.draw[[iter+1]][[1,1]] <- matrix(sapply(1:n, function(i) {
          if (Y_complete[i,] == 1) {
            truncnorm::rtruncnorm(1, a = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          } else {
            truncnorm::rtruncnorm(1, b = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          }
        }), ncol = 1)
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Calculating the joint and individual structure, scaled to the data
  # ---------------------------------------------------------------------------

  # Storing the structures at each Gibbs sampling iteration
  J.draw <- A.draw <- S.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))

  # Calculating the structure for Y at each Gibbs sampling iteration
  EY.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))

  for (iter in 1:nsample) {
    for (s in 1:q) {
      # Calculating the joint structure and scaling by sigma.mat
      J.draw[[iter]][[s,1]] <- (U.train[[iter]][[s,1]] %*% t(V.draw[[iter]][[1,1]])) * sigma.mat[s,1]

      # Calculating the individual structure and scaling by sigma.mat
      A.draw[[iter]][[s,1]] <- (W.train[[iter]][[s,s]] %*% t(Vs.draw[[iter]][[1,s]])) * sigma.mat[s,1]

      # Calculate the overall structure (joint + individual
      S.draw[[iter]][[s,1]] <- J.draw[[iter]][[s,1]] + A.draw[[iter]][[s,1]]
    }

    # Calculate the structure for Y
    if (response_given) {
      if (response_type == "continuous") {
        EY.draw[[iter]][[1,1]] <- VStar.draw[[iter]][[1,1]] %*% beta.train[[iter]][[1,1]]
      }

      if (response_type == "binary") {
        EY.draw[[iter]][[1,1]] <- pnorm(VStar.draw[[iter]][[1,1]] %*% beta.train[[iter]][[1,1]])
      }
    }
  }

  # Return
  list(test_data = test_data, # Returning the test data
       Y_test = Y_test, # Return the test response vector
       J.draw = J.draw, A.draw = A.draw, S.draw = S.draw, EY.draw = EY.draw, # Underlying structure
       V.draw = V.draw, U.train = U.train, W.train = W.train, Vs.draw = Vs.draw, # Components of the structure
       VStar.draw = VStar.draw, # Components that predict Y,
       ranks = c(r, r.vec), # Ranks
       tau2.train = tau2.train, beta.train = beta.train) # Regression parameters

}

#' summarize_factors
#'
#' Calculate posterior summaries of aligned estimated factors from BSFP.
#'
#' @param data A matrix of lists or a list of matrices that share the same number of
#' columns. The matrices must be oriented in pxn orientation. May contain NAs if
#' there are missing values in the dataset.
#' @param Y A matrix of lists or a nx1 matrix of continuous or binary outcome.
#' May be NULL if no outcome is given. May contain NAs if there are missing outcomes.
#' @param iters_burnin (vector): indices for posterior samples after burn-in
#' @param aligned_results (list): results from match_align_bsfp
#' @param ranks Estimated joint and individual ranks from BSFP
#' @param tau2.draw Posterior samples for error variance in Y if response is given
#' @param Xm.draw Imputed values for missing values in X
#' @param Ym.draw Imputed values for unobserved outcomes in Y
#'
#' @details Generate posterior summaries (posterior mean, 95\% credible intervals)
#' for the estimated factors from the joint and individual structures. May then be
#' used to plot the contributions of each biomarker to a factor (loadings) and the
#' expression levels of each factor across samples (scores). This function is
#' intended to be used after alignment using \code{match_align_bsfp()}
#'
#' @return Returns a list of posterior summaries for the joint and individual loadings and scores,
#' regression coefficients, estimated error variance in a continuous response, missing data
#' imputation for missing values in the sources and in the outcome.
#'
#' This list contains the following elements:
#' \item{joint.scores.summary}{List of matrices for each joint factor. Each matrix has \eqn{n} rows
#' for each sample. The first column is the posterior mean score after burn-in for a given sample.
#' Then gives the lower and upper bounds for the 95\% credible interval. If no joint structure
#' was estimated, then returns \code{NULL}.}
#' \item{joint.loadings.summary}{List of lists for each joint factor. Each list contains \eqn{q} lists
#' for the contributions of each biomarker from each source to the given factor. Provides the
#' posterior mean and 95\% credible interval. If no joint structure was estimated, then returns \code{NULL}.}
#' \item{individual.scores.summary}{List of length \eqn{q} with inners lists for each individual factor. Summarizes
#' the score for each sample for each factor. Provides posterior mean and 95\% credible interval. If no individual
#' structure from a source was estimated, returns \code{NULL}. }
#' \item{individual.loadings.summary}{List of length \eqn{q} with inners lists for each individual factor.
#' Summarizes the loadings for each biomarker from each source to the given factor. Provides posterior mean and 95\% credible interval. If no individual
#' structure from a source was estimated, returns \code{NULL}.}
#' \item{joint.betas.summary}{Summarizes the regression coefficients for the joint factors. Provides
#' the posterior mean and 95\% credible interval. If no joint structure was estimated, returns \code{NULL}.}
#' \item{individual.betas.summary}{List of length \eqn{q} with summaries for the regression coefficients
#' of each individual factor from each source. Provides
#' the posterior mean and 95\% credible interval. If no individual structure was estimated for a source,
#' returns \code{NULL}.}
#' \item{Xm.summary}{Returns a list of length \eqn{q} with posterior summaries (mean and 95\% credible interval)
#' for each missing sample from each source. Rownames correspond to the source and the sample index. The index
#' counts from the upper leftmost entry of the matrix, i.e. the index corresponds to \code{which(is.na(data[[s,1]]))}
#' for \eqn{s=1,\dots,q}.}
#' \item{Ym.summary}{Summarizes the imputed values for unobserved outcomes. Provides the posterior mean
#' and 95\% credible interval. Indexes from \eqn{1,\dots, n}. If no unobserved outcomes, returns \code{NULL}.}
#' \item{ranks}{Vector of length \eqn{q+1} of estimated joint and individual ranks. ranks[1] contains the joint rank.
#' ranks[2:(q+1)] contains individual ranks.}
#' \item{tau2.summary}{Posterior mean (mean and 95\% credible interval) of estimated error variance in
#' continuous outcome, if given. If no continuous outcome was used, returns \code{NULL}.}
#' @export
#'
#' @examples
#' # Setting up the data
#' n <- 50
#' p.vec <- c(75, 100)
#' q <- 2
#'
#' # Setting up the model parameters
#' true_params <- list(error_vars = c(1,1),
#'                     joint_var = 1,
#'                    indiv_vars = c(1,1),
#'                    beta_vars = c(1, 1, rep(1, q)),
#'                    response_vars = c(shape = 1, rate = 1))
#'
#' # Choose ranks
#' r <- 3
#' r.vec <- c(3, 3)
#' ranks <- c(r, r.vec)
#'
#' # Number of posterior sampling iterations
#' nsample <- 1000
#' burnin <- nsample/2
#' iters_burnin <- (burnin+1):nsample
#'
#' # Generate data
#' data.c1 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE)
#'
#' # Run BSFP for 1000 iterations
#' bsfp.c1 <- bsfp(data = data.c1$data, Y = data.c1$Y, nsample = nsample)
#'
#' # Run the alignment algorithm
#' alignment.c1 <- match_align_bsfp(BSFP.fit = bsfp.c1, y = data.c1$Y,
#'                                  model_params = bsfp.c1$model_params,
#'                                  p.vec = p.vec, iters_burnin = iters_burnin)
#'
#' # Summarize aligned factors
#' summary.aligned.c1 <- summarize_factors(data = data.c1$data, Y = data.c1$Y,
#'                                         iters_burnin = iters_burnin,
#'                                         aligned_results = alignment.c1,
#'                                         ranks = bsfp.c1$ranks, tau2.draw = bsfp.c1$tau2.draw)

summarize_factors <- function(data, Y = NULL, iters_burnin,
                              aligned_results, ranks, tau2.draw = NULL, Xm.draw = NULL, Ym.draw = NULL) {

  # ---------------------------------------------------------------------------
  # Save the aligned results
  # ---------------------------------------------------------------------------

  joint.scores.final <- aligned_results$joint.scores.final
  joint.loadings.final <- aligned_results$joint.loadings.final
  joint.betas.final <- aligned_results$joint.betas.final

  individual.scores.final <- aligned_results$individual.scores.final
  individual.loadings.final <- aligned_results$individual.loadings.final
  individual.betas.final <- aligned_results$individual.betas.final

  # ---------------------------------------------------------------------------
  # Save data attributes
  # ---------------------------------------------------------------------------

  # Was the data input as a list?
  if (!("matrix" %in% class(data))) {

    # Save the number of sources
    q <- length(data)

    # Initialize new data matrix
    new_data <- matrix(list(), nrow = q, ncol = 1)

    # Add in the sources
    for (s in 1:q) {
      new_data[[s,1]] <- data[[s]]
    }

    # Rename the new version of the data
    data <- new_data
  }

  # If a response was given, was it input as a matrix of lists or as a matrix?
  if (!is.null(Y)) {
    if (class(Y[1,1]) != "list") {

      # Create a new version of Y
      new_Y <- matrix(list(), nrow = 1, ncol = 1)

      # Input the response
      new_Y[[1,1]] <- Y

      # Return to the Y variable
      Y <- new_Y
    }
  }

  # ---------------------------------------------------------------------------
  # Extracting the dimensions
  # ---------------------------------------------------------------------------

  q <- nrow(data) # Number of sources
  p.vec <- apply(data, 1, function(source) nrow(source[[1]])) # Number of features per source
  p <- sum(p.vec) # Total number of features
  n <- ncol(data[[1,1]]) # Number of subjects

  # Initialize the indices of features in each source
  p.ind <- lapply(1:q, function(s) {
    if (s == 1) {
      1:p.vec[s]
    } else {
      (cumsum(p.vec)[s-1] + 1):cumsum(p.vec)[s]
    }
  })

  # Rank indices
  rank.inds <- lapply(1:(q+1), function(s) {

    # For any structure,
    if (ranks[s] == 0) {
      NULL
    }

    # For joint structure
    else if (s == 1) {
      if (ranks[s] > 0) {
        1:ranks[s]
      }
    }

    # For each individual structure
    else if (s > 1) {
      (sum(ranks[1:(s-1)])+1):sum(ranks[1:s])
    }

  })

  # Regression coefficient indices
  beta.ind <- lapply(1:(q+1), function(s) {

    if (is.null(rank.inds[[s]])) {
      NULL
    }

    else {
      rank.inds[[s]] + 1
    }
  })

  # ---------------------------------------------------------------------------
  # Is there a response vector?
  # ---------------------------------------------------------------------------

  response_given <- !is.null(Y[[1,1]])

  # If so, what kind of response is it?
  if (response_given) {
    Y <- matrix(unlist(Y))

    response_type <- if (all(unique(Y) %in% c(0, 1, NA))) "binary" else "continuous"

    # If there is a response, is there missingness in the outcome?
    missingness_in_response <- any(is.na(Y))

    # Which entries are missing?
    missing_obs_Y <- which(is.na(Y))
  }

  # ---------------------------------------------------------------------------
  # Check for missingness in data
  # ---------------------------------------------------------------------------

  # Check for missingness
  missingness_in_data <- any(sapply(data[,1], function(source) any(is.na(source))))

  # Which entries are missing?
  missing_obs <- lapply(data[,1], function(source) which(is.na(source)))

  # ---------------------------------------------------------------------------
  # Joint factors
  # ---------------------------------------------------------------------------

  if (ranks[1] > 0) {
    # Scores
    joint.scores.summary <- lapply(1:ranks[1], function(r) {
      score.matrix <- do.call(cbind, lapply(joint.scores.final, function(iter) iter[,r,drop=FALSE]))
      summary.mat <- cbind.data.frame(Post.Mean = rowMeans(score.matrix),
                                      Lower.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.025)),
                                      Upper.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.975)))
      rownames(summary.mat) <- paste("Sample", 1:n)
      summary.mat
    })
    names(joint.scores.summary) <- paste0("Joint.Factor.", 1:ranks[1])

    # Loadings
    joint.loadings.summary <- lapply(1:q, function(s) {
      lapply(1:ranks[1], function(r) {
        load.matrix <- do.call(cbind, lapply(joint.loadings.final, function(iter) iter[p.ind[[s]],r,drop=FALSE]))
        summary.mat <- cbind.data.frame(Post.Mean = rowMeans(load.matrix),
                                        Lower.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.025)),
                                        Upper.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.975)))
        rownames(summary.mat) <- paste("Source", s, "Feature", 1:p.vec[s])
        summary.mat
      })
    })

    # Regression coefficients
    if (response_given) {
      joint.betas.summary <- lapply(1:ranks[1], function(r) {
        beta.mat <- do.call(cbind, lapply(joint.betas.final, function(iter) iter[r,,drop=FALSE]))
        summary.mat <- cbind.data.frame(Post.Mean = rowMeans(beta.mat),
                                        Lower.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.025)),
                                        Upper.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.975)))
        summary.mat
      })
      joint.betas.summary <- do.call(rbind, joint.betas.summary)
      rownames(joint.betas.summary) <- paste0("Joint.Factor.Regression.Coefficient.", 1:ranks[1])
    }

    if (!response_given) {
      joint.betas.summary <- NULL
    }
  }

  if (ranks[1] == 0) {
    # Scores
    score.matrix <- do.call(cbind, lapply(joint.scores.final, function(iter) iter[,1,drop=FALSE]))
    joint.scores.summary <- cbind.data.frame(Post.Mean = rowMeans(score.matrix),
                                             Lower.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.025)),
                                             Upper.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.975)))
    rownames(joint.scores.summary) <- paste("Sample", 1:n)
    names(joint.scores.summary) <- paste0("Joint.Factor.", 1:ranks[1])

    # Loadings
    joint.loadings.summary <- lapply(1:q, function(s) {
      load.matrix <- do.call(cbind, lapply(joint.loadings.final, function(iter) iter[p.ind[[s]],1,drop=FALSE]))
      summary.mat <- cbind.data.frame(Post.Mean = rowMeans(load.matrix),
                                      Lower.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.025)),
                                      Upper.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.975)))
      rownames(summary.mat) <- paste("Source", s, "Feature", 1:p.vec[s])
      summary.mat
    })

    # Regression coefficients
    if (response_given) {
      beta.mat <- do.call(cbind, lapply(joint.betas.final, function(iter) iter[1,,drop=FALSE]))
      joint.betas.summary <- cbind.data.frame(Post.Mean = rowMeans(beta.mat),
                                              Lower.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.025)),
                                              Upper.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.975)))
      joint.betas.summary
    }

    if (!response_given) {
      joint.betas.summary <- NULL
    }

  }

  # ---------------------------------------------------------------------------
  # Individual factors
  # ---------------------------------------------------------------------------

  # Scores
  individual.scores.summary <- lapply(1:q, function(s) {

    if (ranks[s+1] > 0) {
      source.scores.list <- lapply(1:ranks[s+1], function(rs) {
        score.matrix <- do.call(cbind, lapply(individual.scores.final[[s]], function(iter) iter[,rs,drop=FALSE]))
        summary.mat <- cbind.data.frame(Post.Mean = rowMeans(score.matrix),
                                        Lower.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.025)),
                                        Upper.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.975)))
        rownames(summary.mat) <- paste("Sample", 1:n)
        summary.mat
      })
      names(source.scores.list) <- paste0("Source.", s, ".Individual.Factor.", 1:ranks[s+1])
      source.scores.list
    } else {
      score.matrix <- do.call(cbind, lapply(individual.scores.final[[s]], function(iter) iter[,1,drop=FALSE]))
      summary.mat <- cbind.data.frame(Post.Mean = rowMeans(score.matrix),
                                      Lower.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.025)),
                                      Upper.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.975)))
      rownames(summary.mat) <- paste("Sample", 1:n)
      summary.mat
    }

  })
  names(individual.scores.summary) <- paste0("Source.", 1:q)

  # Loadings
  individual.loadings.summary <- lapply(1:q, function(s) {

    if (ranks[s+1] > 0) {
      source.load.list <- lapply(1:ranks[s+1], function(rs) {
        load.matrix <- do.call(cbind, lapply(individual.loadings.final[[s]], function(iter) iter[,rs,drop=FALSE]))
        summary.mat <- cbind.data.frame(Post.Mean = rowMeans(load.matrix),
                                        Lower.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.025)),
                                        Upper.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.975)))
        rownames(summary.mat) <- paste("Source", s, "Feature", 1:p.vec[s])
        summary.mat
      })
      names(source.load.list) <- paste0("Source.", s, ".Individual.Factor.", 1:ranks[s+1])
      source.load.list
    } else {
      load.matrix <- do.call(cbind, lapply(individual.loadings.final[[s]], function(iter) iter[,1,drop=FALSE]))
      summary.mat <- cbind.data.frame(Post.Mean = rowMeans(load.matrix),
                                      Lower.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.025)),
                                      Upper.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.975)))
      rownames(summary.mat) <- paste("Source", s, "Feature", 1:p.vec[s])
      summary.mat
    }

  })
  names(individual.loadings.summary) <- paste0("Source.", 1:q)

  # Regression coefficients
  if (response_given) {
    individual.betas.summary <- lapply(1:q, function(s) {

      if (ranks[s+1] > 0) {
        source.beta.list <- lapply(1:ranks[s+1], function(rs) {
          beta.mat <- do.call(cbind, lapply(individual.betas.final[[s]], function(iter) iter[rs,,drop=FALSE]))
          summary.mat <- cbind.data.frame(Post.Mean = rowMeans(beta.mat),
                                          Lower.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.025)),
                                          Upper.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.975)))
          summary.mat
        })
        source.beta.list <- do.call(rbind, source.beta.list)
        rownames(source.beta.list) <- paste0("Source.", s, ".Individual.Factor.Regression.Coefficient.", 1:ranks[s+1])
        source.beta.list
      } else {
        beta.mat <- do.call(cbind, lapply(individual.betas.final[[s]], function(iter) iter[1,,drop=FALSE]))
        summary.mat <- cbind.data.frame(Post.Mean = rowMeans(beta.mat),
                                        Lower.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.025)),
                                        Upper.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.975)))
        summary.mat
      }

    })
    names(individual.betas.summary) <- paste0("Source.", 1:q)
  }

  if (!response_given) {
    individual.betas.summary <- NULL
  }

  # ---------------------------------------------------------------------------
  # Summarizing the outcome variance
  # ---------------------------------------------------------------------------

  if (response_given) {
    if (!is.null(tau2.draw)) {
      tau2.summary <- cbind.data.frame(Post.Mean = mean(unlist(tau2.draw)[iters_burnin]),
                                       Lower.95.CI = quantile(unlist(tau2.draw)[iters_burnin], 0.025),
                                       Upper.95.CI = quantile(unlist(tau2.draw)[iters_burnin], 0.975))
    }
  }

  if (!response_given | is.null(tau2.draw)) {
    tau2.summary <- NULL
  }

  # ---------------------------------------------------------------------------
  # Summarizing the imputed values
  # ---------------------------------------------------------------------------

  # For each source
  Xm.summary <- lapply(1:q, function(s) {

    if (length(missing_obs[[s]]) > 0) {
      # Create matrix of imputed values for each missing sample
      Xm.s <- do.call(cbind, lapply(iters_burnin, function(iter) {
        Xm.draw[[iter]][[s,1]]
      }))

      # Summarize the posterior mean and 95% CI
      Xm.s.summary <- cbind.data.frame(Post.Mean = rowMeans(Xm.s),
                                       Lower.95.CI = apply(Xm.s, 1, function(row) quantile(row, 0.025)),
                                       Upper.95.CI = apply(Xm.s, 1, function(row) quantile(row, 0.975)))

      # Add rownames
      rownames(Xm.s.summary) <- paste0("Source ", s, " Sample ", missing_obs[[s]])
      Xm.s.summary
    } else {
      NULL
    }
  })

  # For the response
  if (response_given) {
    if (missingness_in_response) {
      Ym.mat <- do.call(cbind, lapply(iters_burnin, function(iter) Ym.draw[[iter]][[1,1]]))
      Ym.summary <- cbind.data.frame(Post.Mean = rowMeans(Ym.mat),
                                     Lower.95.CI = apply(Ym.mat, 1, function(row) quantile(row, 0.025)),
                                     Upper.95.CI = apply(Ym.mat, 1, function(row) quantile(row, 0.975)))
      rownames(Ym.summary) <- paste0("Sample ", missing_obs_Y)
    }

    if (!missingness_in_response) {
      Ym.summary <- NULL
    }
  }

  if (!response_given) {
    Ym.summary <- NULL
  }


  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------

  list(
    # Summaries of the joint factors
    joint.scores.summary = joint.scores.summary,
    joint.loadings.summary = joint.loadings.summary,

    # Summaries of the individual factors
    individual.scores.summary = individual.scores.summary,
    individual.loadings.summary = individual.loadings.summary,

    # Summaries of the factor contributions in predictive model
    joint.betas.summary = joint.betas.summary,
    individual.betas.summary = individual.betas.summary,

    # Summaries of the imputed values
    Xm.summary = Xm.summary, Ym.summary = Ym.summary, # Missing data imputation

    # Miscellaneous
    ranks = ranks, # Ranks
    tau2.summary = tau2.summary) # Regression parameters

}

#' plot_summaries
#'
#' Plot the posterior summaries (posterior mean, 95% credible interval) for the estimated
#' scores, loadings, and regression coefficients.
#'
#' @param summaries Output from \code{summarize_factors()} function
#' @param structure (string) specify one of "joint" or "individual"
#' @param output (string) specify one of "scores", "loadings", "betas"
#' @param source (vector of ints) specify source by index from 1:q or NULL if plotting the joint scores
#' @param source.name (string) Specify the name of the source, e.g. "Expression" or "Metabolome"
#' @param sample.labels (vector of strings) NULL if \code{output != scores}. Otherwise, labels for plotting sample scores, e.g. c(1, "Control", "Control", ...). Could be the outcome used in BSFP. If NULL, plot will highlight credible intervals which don't contain 0.
#' @param biomarker.names (vector of strings) NULL if \code{structure != "individual}. Otherwise, labels for plotting the biomarker loadings, e.g. c("Gene 1", "Gene 2", ...). If NULL, plot will highlight credible intervals which don't contain 0.
#' @param label.x (Boolean) Should the x-axis labels be plotted? Default is FALSE. There may be many labels which makes the plot messy.
#'
#' @details Plotting function for the posterior means and credible intervals generated from \code{summarize_factors}.
#' Automatically orders the scores/loadings/regression coefficients by posterior mean for easy visualization.
#' User can specify loadings or scores to be colored based on a factor using \code{biomarker.labels} or \code{sample.labels}, respectively.
#' Outputs a list of plots for scores and loadings. See example for how to conveniently view each plot.
#'
#' @return Returns a list of plots for each factor if plotting loadings or scores. If plotting regression coefficients,
#' outputs a list with 1 plot only.
#'
#' @export
#'
#' @examples
#' #' # Setting up the data
#' n <- 50
#' p.vec <- c(75, 100)
#' q <- 2
#'
#' # Setting up the model parameters
#' true_params <- list(error_vars = c(1,1),
#'                     joint_var = 1,
#'                    indiv_vars = c(1,1),
#'                    beta_vars = c(1, 1, rep(1, q)),
#'                    response_vars = c(shape = 1, rate = 1))
#'
#' # Choose ranks
#' r <- 3
#' r.vec <- c(3, 3)
#' ranks <- c(r, r.vec)
#'
#' # Number of posterior sampling iterations
#' nsample <- 1000
#' burnin <- nsample/2
#' iters_burnin <- (burnin+1):nsample
#'
#' # Generate data
#' data.c1 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE)
#'
#' # Run BSFP for 1000 iterations
#' bsfp.c1 <- bsfp(data = data.c1$data, Y = data.c1$Y, nsample = nsample)
#'
#' # Run the alignment algorithm
#' alignment.c1 <- match_align_bsfp(BSFP.fit = bsfp.c1, y = data.c1$Y,
#'                                  model_params = bsfp.c1$model_params,
#'                                  p.vec = p.vec, iters_burnin = iters_burnin)
#'
#' # Summarize aligned factors
#' summary.aligned.c1 <- summarize_factors(data = data.c1$data, Y = data.c1$Y,
#'                                         iters_burnin = iters_burnin,
#'                                         aligned_results = alignment.c1,
#'                                         ranks = bsfp.c1$ranks, tau2.draw = bsfp.c1$tau2.draw)
#' # Plot summaries
#' plots.joint.scores <- plot_summaries(summary.aligned.c1, structure = "joint", output = "scores")
#' plots.joint.loadings.source1 <- plot_summaries(summary.aligned.c1, structure = "joint", output = "loadings", source = 1, xlab.name = "Genes")
#' plots.joint.betas <- plot_summaries(summary.aligned.c1, structure = "joint", output = "betas")
#'
#' plots.individual.scores.source2 <- plot_summaries(summary.aligned.c1, structure = "individual", output = "scores", source = 2)
#' plots.individual.loadings.source2 <- plot_summaries(summary.aligned.c1, structure = "individual", output = "loadings", source = 2, xlab.name = "Proteins")
#' plots.individuaul.betas.source2 <- plot_summaries(summary.aligned.c1, structure = "individual", output = "betas", source = 2)
#'
#' # View one at a time
#' plots.joint.scores[[1]] # Joint factor 1
#' plots.joint.loadings.source1[[2]] # Joint factor 2
#' plots.joint.betas[[1]] # All regression coefficients for joint factors
#'
#' plots.individual.scores.source2[[1]] # Scores for individual factor 1 from source 2
#' plots.individual.loadings.source2[[2]] # Loadings for individual factor 2 from source 2
#' plots.individuaul.betas.source2[[1]] # All regression coefficients for individual factors for source 2
#'
#' # OR, output to a pdf (will output to current working directory)
#' pdf("Joint_Scores_BSFP.pdf")
#' plots.joint.scores
#' dev.off()
plot_summaries <- function(summaries, structure, output, source = NULL, xlab.name = NULL,
                           sample.labels = NULL, biomarker.labels = NULL, label.x = FALSE) {

  # summaries = output from summarize_factors
  # structure = "joint" or "individual"
  # output = "scores", "loadings", "betas"
  # source = 1:q or NULL if plotting joint scores
  # biomarker.names = c("Source 1 name", "Source 2 name", ...) e.g. c("Gene", "Metabolite", ...) or NULL if structure != "individual"
  # grouping = label for plotting scores
  # label.x = Boolean. Should x-axis be labeled? Could be messy with many samples/biomarkers/factors
  # sample.labels = labels for samples if desired. NULL if output != scores or if not wanted. Default will plot CIs not containing 0.
  # biomarker.labels = labels for biomarkers if desired. NULL if output != loadings or if not wanted. Default will plot CIs not containing 0.

  # Paste together structure and output names
  structure.output <- paste0(structure, ".", output)

  # Check if source was specified if not plotting joint scores
  if (structure.output != "joint.scores" & structure.output != "joint.betas" & is.null(source)) stop("You must specify the source argument (choose from 1:q where q is the number of sources)")


  # Check if an x-axis label was specified
  if (output == "loadings" & is.null(xlab.name)) stop("You must specify a name for the x-axis label using xlab.name, e.g. xlab.name = 'Gene Loadings' ")


  # Select the object from the summary
  obj.summary <- summaries[grepl(structure.output, names(summaries))][[1]]

  # If structure == "individual", save the corresponding source
  if (structure == "individual" | output == "loadings") {
    obj.summary <- obj.summary[[source]]
  }

  # Initialize a list of plots
  if (output != "betas") {
    plot.list <- lapply(1:length(obj.summary), function(i) list())
  }

  if (output == "betas") {
    obj.summary <- list(obj.summary)
    plot.list <- list(1)
  }

  # Iterate through the factors and create plots
  for (i in 1:length(plot.list)) {

    # Add the sample/biomarker ID
    obj.summary.i <- cbind.data.frame(ID = factor(rownames(obj.summary[[i]])),
                                      obj.summary[[i]])

    # Add labels for the bimoarkers if given
    if (!is.null(biomarker.labels)) {
      obj.summary.i <- cbind.data.frame(ID = factor(rownames(obj.summary[[i]])),
                                        label = factor(biomarker.labels),
                                        obj.summary[[i]])
    }

    # If labels for the samples are given
    if (!is.null(sample.labels)) {
      obj.summary.i <- cbind.data.frame(ID = factor(rownames(obj.summary[[i]])),
                                        label = factor(sample.labels),
                                        obj.summary[[i]])
    }

    # Order by mean
    obj.summary.i <- obj.summary.i[order(obj.summary.i$Post.Mean, decreasing = TRUE),]

    # Add whether or not CI contains 0
    if (is.null(biomarker.labels) & is.null(sample.labels)) {
      obj.summary.i.sig <- obj.summary.i %>%
        filter(sign(Lower.95.CI) == sign(Upper.95.CI))
    }

    # Color the CIs by the non-reference label if given
    if (!is.null(biomarker.labels)) {
      obj.summary.i.sig <- obj.summary.i %>%
        filter(label == levels(obj.summary.i$label)[length(levels(obj.summary.i$label))])
    }

    if (!is.null(sample.labels)) {
      obj.summary.i.sig <- obj.summary.i %>%
        filter(label == levels(obj.summary.i$label)[length(levels(obj.summary.i$label))])
    }

    # Set the title
    structure.title <- ifelse(structure == "joint", "Joint", "Individual")
    output.title <- "Scores"
    output.title <- ifelse(output != "scores", ifelse(output == "loadings", "Loadings", "Regression Coefficients"), "Scores")

    # Depending on the output type, change the title accordingly
    # If output is not the regression coefficients, state specific factor
    if (output != "betas") {
      title <- paste0(structure.title, " Factor ", i, " - ", output.title, " and Credible Intervals")
    }

    # If output is the regression coefficients, mention specific structure
    if (output == "betas") {
      title <- paste0(structure.title, " Structure - ", output.title, " and Credible Intervals")
    }

    # Set the x- and y-axis labels
    xlab <- "Sample"
    xlab <- ifelse(output != "scores", ifelse(output == "loadings", xlab.name, "Regression Coefficients"), "Sample")
    ylab <- paste0("Posterior Mean")

    # Plot the results for the scores
    pp <- obj.summary.i %>%
      ggplot(aes(ID, Post.Mean)) +
      geom_point() +
      geom_errorbar(aes(ymin = Lower.95.CI, ymax = Upper.95.CI)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      aes(x = fct_inorder(ID)) +
      geom_errorbar(aes(ymin = Lower.95.CI, ymax = Upper.95.CI),
                    col = "black") +
      geom_point(data = obj.summary.i.sig,
                 aes(ID, Post.Mean), color = "orange") +
      geom_errorbar(data = obj.summary.i.sig,
                    aes(ymin = Lower.95.CI, ymax = Upper.95.CI),
                    col = "orange") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggtitle(title) +
      xlab(xlab) +
      ylab(ylab)

    # If x-axis labels should be removed
    if (!label.x) {
      pp <- pp + theme(axis.text.x=element_blank())
    }

    # Add figure to list
    plot.list[[i]] <- pp
  }

  # Return
  plot.list
}

#' Log-joint density of Bayesian Simultaneous Factorization and Prediction
#'
#' Calculate the log-joint density of the estimated model at a given posterior
#' sampling iteration to facilitate checking convergence.
#'
#' @param data A matrix of lists or a list of matrices that share the same number of
#' columns. The matrices must be oriented in pxn orientation. May contain NAs if
#' there are missing values in the dataset.
#' @param Y A matrix of lists or a nx1 matrix of continuous or binary outcome.
#' May be NULL if no outcome is given. May contain NAs if there are missing outcomes.
#' @param U.iter Posterior draw for joint loadings, U, at a given iteration
#' @param V.iter Posterior draw for joint scores, V, at a given iteration
#' @param W.iter Posterior draw for individual loadings, W, at a given iteration
#' @param Vs.iter Posterior draw for individual scores, Vs, at a given iteration
#' @param model_params Model parameters used in BSFP fit. Usually given in BSFP output.
#' @param ranks Estimated joint and individual ranks from BSFP
#' @param beta.iter Posterior draw for regression coefficients at a given iteration
#' @param tau2.iter Posterior draw for estimated error variance in Y at a given iteration
#' @param Xm.iter Imputed values for unobserved X values at a given iteration
#' @param Ym.iter Imputed values for unobserved Y values at a given iteration
#'
#' @details Calculates the log-joint density of the BSFP model at each posterior sampling iteration.
#' Use to determine if model appears to have converged. We expect to see the log-joint density
#' increase and then plateau after convergence. Upon convergence, the log-joint density should
#' bounce around randomly.
#'
#' @returns Returns a single value representing the log-joint density at a given posteiror
#' sampling iteration.
#'
#' @export
#'
#' @examples
#' # Setting up the data
#' n <- 50
#' p.vec <- c(75, 100)
#' q <- 2
#'
#' # Setting up the model parameters
#' true_params <- list(error_vars = c(1,1),
#'                     joint_var = 1,
#'                     indiv_vars = c(1,1),
#'                     beta_vars = c(1, 1, rep(1, q)),
#'                     response_vars = c(shape = 1, rate = 1))
#'
#' # Choose ranks
#' r <- 3
#' r.vec <- c(3, 3)
#' ranks <- c(r, r.vec)
#'
#' # Number of posterior sampling iterations
#' nsample <- 1000
#' burnin <- nsample/2
#' iters_burnin <- (burnin+1):nsample
#'
#' # Generate data
#' data.c1 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE)
#'
#' # Run BSFP for 1000 iterations
#' bsfp.c1 <- bsfp(data = data.c1$data, Y = data.c1$Y, nsample = nsample)
#'
#' # Check convergence
#' log_joint_density_by_iter.c1 <- sapply(1:nsample, function(iter) {
#'   log_joint_density(data = data.c1$data, Y = data.c1$Y,
#'                     U.iter = bsfp.c1$U.draw[[iter]],
#'                     V.iter = bsfp.c1$V.draw[[iter]],
#'                     W.iter = bsfp.c1$W.draw[[iter]],
#'                     Vs.iter = bsfp.c1$Vs.draw[[iter]],
#'                     model_params = bsfp.c1$model_params,
#'                    ranks = bsfp.c1$ranks,
#'                    beta.iter = bsfp.c1$beta.draw[[iter]],
#'                    tau2.iter = bsfp.c1$tau2.draw[[iter]])
#' })

log_joint_density <- function(data, Y = NULL, U.iter, V.iter, W.iter, Vs.iter, model_params, ranks, beta.iter = NULL, tau2.iter = NULL, Xm.iter = NULL, Ym.iter = NULL) {

  # How many sources are there?
  q <- nrow(data)

  # How many observations?
  n <- ncol(data[[1,1]])

  # Saving the model parameters
  error_vars <- model_params$error_vars # Error variances
  sigma2_joint <- joint_var <- model_params$joint_var # Variance of joint structure
  sigma2_indiv <- indiv_vars <- model_params$indiv_vars # Variances of individual structure
  beta_vars <- model_params$beta_vars # Variances on betas
  response_vars <- model_params$response_vars; shape <- response_vars[1]; rate <- response_vars[2] # Hyperparameters of variance of response

  r <- ranks[1]
  r.vec <- ranks[-1]
  r_total <- r + sum(r.vec)
  n_beta <- 1 + r_total

  Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
  diag(Sigma_beta) <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))

  # Check if there is a response
  response_given <- !is.null(Y)

  # If there is a response, what type of response is it?
  if (response_given) {
    response_type <- if (all(unique(Y) %in% c(0, 1, NA))) "binary" else "continuous"
  }

  # Check if there is missingness
  missingness_in_data <- any(sapply(data[,1], function(source) any(is.na(source))))

  # Which entries are missing?
  missing_obs <- lapply(data[,1], function(source) which(is.na(source)))

  # If there is missingness, fill in the missing values with the current imputed values
  if (missingness_in_data) {
    # Creating the completed matrices.
    X_complete <- data

    # Fill in the completed matrices with the imputed values
    for (s in 1:q) {
      X_complete[[s,1]][missing_obs[[s]]] <- Xm.iter[[s,1]]
    }
  }

  # If there is no missingness, rename the data
  if (!missingness_in_data) {
    X_complete <- data
  }

  # Check if there is missingness in the response
  missingness_in_response <- any(is.na(Y[[1,1]]))

  # Which entries are missing?
  missing_obs_Y <- which(is.na(Y[[1,1]]))

  # If there is missingness in the response, fill in the missing values with the current imputed values
  if (missingness_in_response) {
    # Create a completed response
    Y_complete <- Y

    # Fill in the missing values
    Y_complete[[1,1]][missing_obs_Y] <- Ym.iter[[1,1]]
  }

  # If there is no missingness, rename the response
  if (!missingness_in_response) {
    Y_complete <- Y
  }

  # Saving the contributions of each term to the joint density
  like <- 0

  # Contribution of V to the joint density
  like <- like + sum(sapply(1:r, function(rs) {
    dnorm(V.iter[[1,1]][,rs], mean = 0, sd = sqrt(sigma2_joint), log = TRUE)
  }))

  for (s in 1:q) {
    # Contribution of the observed data to the joint density
    data_s <- X_complete[[s,1]]
    like <- like + sum(sapply(1:n, function(i) {
      dnorm(data_s[,i], mean = (U.iter[[s,1]] %*% t(V.iter[[1,1]]) + W.iter[[s,s]] %*% t(Vs.iter[[1,s]]))[,i], sd = sqrt(error_vars[s]), log = TRUE)
    }))

    # Contribution of Us to the joint density
    like <- like + sum(sapply(1:r, function(rs) {
      dnorm(U.iter[[s,1]][,rs], mean = 0, sd = sqrt(sigma2_joint), log = TRUE)
    }))

    # Contribution of Ws to the joint density
    like <- like + sum(sapply(1:r.vec[s], function(rs) {
      dnorm(W.iter[[s,s]][,rs], mean = 0, sd = sqrt(sigma2_indiv[s]), log = TRUE)
    }))

    # Contribution of Vs to the joint density
    like <- like + sum(sapply(1:r.vec[s], function(rs) {
      dnorm(Vs.iter[[1,s]][,rs], mean = 0, sd = sqrt(sigma2_indiv[s]), log = TRUE)
    }))

    # If there is a response
    if (response_given) {
      VStar.iter <- cbind(1, do.call(cbind, V.iter), do.call(cbind, Vs.iter))

      # The contribution of beta to the joint density
      like <- like + sum(sapply(1:n_beta, function(rs) {
        dnorm(beta.iter[[1,1]][rs,], mean = 0, sd = sqrt(Sigma_beta[rs,rs]), log = TRUE)
      }))

      if (response_type == "continuous") {
        # The contribution of the observed response to the joint density
        like <- like + sum(sapply(1:n, function(i) {
          dnorm(Y_complete[[1,1]][i,], mean = (VStar.iter %*% beta.iter[[1,1]])[i,], sd = sqrt(tau2.iter[[1,1]]), log = TRUE)
        }))

        # The contribution of tau2 to the joint density
        like <- like + log(invgamma::dinvgamma(tau2.iter[[1,1]], shape = shape, scale = 1/rate))
      }

      if (response_type == "binary") {
        # The contribution of the observed response to the joint density
        like <- like + sum(sapply(1:n, function(i) {
          dbinom(Y_complete[[1,1]][i,], size = 1, prob = pnorm(VStar.iter %*% beta.iter[[1,1]])[i,], log = TRUE)
        }))
      }
    }
  }

  # Return
  like
}

#' bsfp_data
#'
#' Generate simulated data according to BSFP model for testing and simulation purposes.
#'
#' @param p.vec (vector) vector of with number of features per source
#' @param n (int) number of samples
#' @param ranks (vector) vector with joint and individual ranks (first index
#' is joint, remaining are individual for each source)
#' @param s2nX (dbl) signal-to-noise ratio in X (ratio of variance in X to variance
#' in error). >1 means higher signal in X, <1 means lower signal in X
#' @param s2nY (dbl) signal-to-noise ratio in Y
#' @param response NULL if no response is desired, "continuous", or "binary"
#' @param missingness NULL if no missingness, or "missingness_in_data" or "missingness_in_response"
#' @param missing_data_type NULL if missingness = NULL; otherwise, = entrywise
#' if randomly sample across sources; columnwise if randomly set samples to missing;
#'  MNAR is randomly set features below a threshold to missing
#' @param prop_missing NULL if missingness = NULL; otherwise the proportion of
#' entries or columns set to missing
#' @param sparsity (Boolean) TRUE if generate regression coefficients under
#' spike-and-slab prior, FALSE otherwise
#' @param identically_zero (Boolean) TRUE if generating response with sparsity
#' and want spike coefficients to be exactly 0
#' @param num_in_spike (vector) number of coefficients to be generated from spike
#' if generating response with sparsity. NULL if sparsity=FALSE. Should have an
#' integer for joint and each individual structure. May be set to NULL and
#' the number of factors in the spike will be determined randomly.
#'
#' @details Generate data according to the BSFP models. Use for examples or
#' simulations to assess model performance. Can scale data sources and reponse
#' vector to accommodate different levels of signal-to-noise ratios. Can also
#' accommodate different kinds of missingness mechanisms (missing-at-random
#' and missingness not at random). This includes removing entire samples from a
#' source. Sparsity arguments are for testing purposes only and the BSFP model
#' does not accommodate sparsity in the predictive model.
#'
#' @returns Returns a list of generated datasets and true parameter values.
#' \item{data}{Returns the observed data without missingness. Data is in matrix-list format,
#' where each dataset is contained within a list entry in a larger matrix object.}
#' \item{Y}{Returns the observed response vector. This will not contain any missingness.}
#' \item{missing_data}{Returns a version of the observed data with missingness induced.
#' Missing values are represented with NAs. If no missingness in the data, this value is \code{NULL}.}
#' \item{missing_obs}{Indices corresponding to missing entries in each dataset. The indices
#' refer to results from \code{which(is.na(data[[s,1]]))} for \eqn{s=1,\dots, q}.
#' If no missingness in the data, this value is \code{NULL}.}
#' \item{Y_missing}{The response vector with missingness induced. Missing values are represented
#' by NAs. If no missingness in the response, this value is \code{NULL}.}
#' \item{s2nX}{User-provided signal-to-noise ratio in the data. If not provided, is \code{NULL}.}
#' \item{s2nX_coef}{Coefficient used to scale the data structure to have the given \code{s2nX}.
#' If no \code{s2nX} is specified, this value is \code{NULL}.}
#' \item{s2nY}{User-provided signal-to-noise ratio in the response If not provided, is \code{NULL}.}
#' \item{s2nY_coef}{Coefficient used to scale \eqn{\mathbb{E}(y|X)} to have the given \code{s2nX}.
#' If no \code{s2nY} is specified, this value is \code{NULL}.}
#' \item{joint.structure}{Returns the true underlying joint structure for each source in matrix-list form.
#' The joint structure is identifiable by BSFP. }
#' \item{indiv.structure}{Returns the true underlying individual structure for each source in matrix-list form.
#' The individual structure is identifiable by BSFP.}
#' \item{V}{The joint scores used to construct the joint structure. This is not identifiable by
#' BSFP. Access the scores by \code{V[[1,1]]}}
#' \item{U}{The joint loadings used to construct the joint structure. This is not identifiable by
#' BSFP. Access the loadings for source \eqn{s} by \code{U[[s,1]]}.}
#' \item{Vs}{The individual scores used to construct the individual structure. This is not identifiable
#' by BSFP. Access the scores for source \eqn{s} by \code{Vs[[1,s]]}.}
#' \item{W}{The individual loadings used to construct the individual structure. This is not identifiable
#' by BSFP. Access the loadings for source \eqn{s} by \code{W[[s,s]]}. \code{W[[s,ss]]} for \eqn{ss \neq s}
#' is equal to 0.}
#' \item{beta}{The intercept and regression coefficients for the joint and individual factors.}
#' \item{EY}{The conditional expectation of the response, i.e. \eqn{\mathbb{E}(y|X)}.}
#' \item{tau2}{The true variance of a continuous response vector. If \code{response="binary"}, this
#' value is \code{NULL}.}
#' \item{gamma}{The true spike-and-slab indicators if response vector is generated with sparsity.
#' 1 reflects inclusion, 0 reflects exclusion.}
#' \item{p.prior}{The prior inclusion probability for the spike-and-slab prior if the response
#' vector is generated with sparsity.}
#'
#' @export
#'
#' @examples
#' # Setting up the data
#' n <- 50
#' p.vec <- c(75, 100)
#' q <- 2
#'
#' # Choose ranks
#' r <- 3
#' r.vec <- c(3, 3)
#' ranks <- c(r, r.vec)
#'
#' # Setting up the model parameters
#' true_params <- list(error_vars = c(1,1),
#'                   joint_var = 1,
#'                     indiv_vars = c(1,1),
#'                     beta_vars = c(1, 1, rep(1, q)),
#'                     response_vars = c(shape = 1, rate = 1))
#'
#' # Generate data with a continuous response and no missingness
#' data.c1 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE)
#'
#' # Generate data with entrywise missingness
#' data.c2 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE,
#' missingness = "missingness_in_data", missing_data_type = "entrywise", prop_missing = 0.1)

bsfp_data <- function(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = NULL, missingness = NULL, missing_data_type = NULL, prop_missing = NULL, sparsity, identically_zero = FALSE, num_in_spike = NULL) {

  # -------------------------------------------------------------------------
  # Setting the dimensions and latent components
  # -------------------------------------------------------------------------

  q <- length(p.vec)
  n_beta <- 1 + sum(ranks)
  r <- ranks[1]
  r.vec <- ranks[-1]

  # Setting the true parameters
  error_vars <- true_params$error_vars
  sigma2_joint <- joint_var <- true_params$joint_var
  sigma2_indiv <- indiv_vars <- true_params$indiv_vars
  beta_vars <- true_params$beta_vars
  response_vars <- true_params$response_vars

  if (!is.null(response_vars)) {
    shape <- response_vars[1]; rate <- response_vars[2]
  }

  # -------------------------------------------------------------------------
  # Generating the underlying structure
  # -------------------------------------------------------------------------

  joint.structure <- indiv.structure <- overall.structure <- matrix(list(), ncol = 1, nrow = q)

  V <- matrix(list(), nrow = 1, ncol = 1)

  if (r > 0) V[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
  if (r == 0) V[[1,1]] <- matrix(0, nrow = n, ncol = 1)

  U <- matrix(list(), nrow = q, ncol = 1)
  Vs <- matrix(list(), nrow = 1, ncol = q)
  W <- matrix(list(), nrow = q, ncol = q)

  E <- matrix(list(), nrow = q, ncol = 1)

  for (s in 1:q) {
    if (r > 0) U[[s,1]] <- matrix(rnorm(p.vec[s]*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = p.vec[s], ncol = r)
    if (r == 0) U[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)

    if (r.vec[s] > 0) {
      Vs[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = n, ncol = r.vec[s])
      W[[s,s]] <- matrix(rnorm(p.vec[s]*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = p.vec[s], ncol = r.vec[s])

      for (ss in 1:q) {
        if (ss != s) {
          if (r.vec[ss] > 0) W[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])

          if (r.vec[ss] == 0) W[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
        }
      }
    }
    if (r.vec[s] == 0) {
      Vs[[1,s]] <- matrix(0, nrow = n, ncol = 1)
      W[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)

      for (ss in 1:q) {
        if (ss != s) {
          if (r.vec[ss] > 0) W[[s,ss]] <- matrix(0, nrow = p.vec[s], ncol = r.vec[ss])

          if (r.vec[ss] == 0) W[[s,ss]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }
      }
    }

    E[[s,1]] <- matrix(rnorm(p.vec[s]*n, sd = sqrt(error_vars[s])), nrow = p.vec[s], ncol = n)

    joint.structure[[s,1]] <- U[[s,1]] %*% t(V[[1,1]])
    indiv.structure[[s,1]] <- W[[s,s]] %*% t(Vs[[1,s]])
    overall.structure[[s,1]] <- joint.structure[[s,1]] + indiv.structure[[s,1]]
  }

  # -------------------------------------------------------------------------
  # Standardizing the variance of the signal in the data
  # -------------------------------------------------------------------------

  if (is.null(s2nX)) {
    s2nX_coef <- NULL
  }

  if (!is.null(s2nX)) {
    # Calculating the scaling coefficient so that the variance of the underlying structure = s2nX * noise variance
    s2nX_coef <- rep(0, q)

    for (s in 1:q) {
      s2nX_coef[s] <- sqrt(s2nX) * sd(c(E[[s,1]]))/sd(c(joint.structure[[s,1]] + indiv.structure[[s,1]]))

      joint.structure[[s,1]] <- s2nX_coef[s] * joint.structure[[s,1]]
      indiv.structure[[s,1]] <- s2nX_coef[s] * indiv.structure[[s,1]]
      overall.structure[[s,1]] <- joint.structure[[s,1]] + indiv.structure[[s,1]]
    }
  }

  # -------------------------------------------------------------------------
  # Calculating the observed data
  # -------------------------------------------------------------------------

  data <- matrix(list(), nrow = q, ncol = 1)

  for (s in 1:q) {
    data[[s,1]] <- joint.structure[[s,1]] + indiv.structure[[s,1]] + E[[s,1]]
  }

  # -------------------------------------------------------------------------
  # Adding a response if desired
  # -------------------------------------------------------------------------

  if (is.null(response)) {
    Y <- EY <- Y_missing <- beta <- tau2 <- gamma <- p.prior <- matrix(list(), nrow = 1, ncol = 1)
    s2nY_coef <- NULL
  }

  if (!is.null(response)) {

    Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
    diag(Sigma_beta) <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))

    if (!sparsity) {
      # Generate betas
      beta <- matrix(list(), nrow = 1, ncol = 1)
      beta[[1,1]] <- matrix(MASS::mvrnorm(1, mu = rep(0, n_beta), Sigma = Sigma_beta), ncol = 1)
      p.prior <- gamma <- matrix(list(), nrow = 1, ncol = 1)
    }

    if (sparsity) {

      p.prior <- matrix(rbeta(1, 1, 1), ncol = 1)

      # If no prior specification on number of coefficients in spike
      if (is.null(num_in_spike)) {
        gamma <- matrix(rbinom(n_beta, size = 1, prob = p.prior), ncol = 1)
        gamma[1,] <- 1 # Always include the intercept
      }

      # If user specifies number of coefficients in spike, randomly choose which
      if (!is.null(num_in_spike)) {
        gamma <- matrix(1, nrow = n_beta, ncol = 1)

        # Randomly choose which components to be in the spike from each joint and individual structure
        spike_inds <- c()

        # Randomly choose joint components for spike
        joint_spike_inds <- sort(sample(x = c(2:r), size = num_in_spike[1], replace = FALSE))
        spike_inds <- c(spike_inds, joint_spike_inds)

        # Randomly choose which individual components for spike
        for (s in 1:q) {
          if (s == 1) {
            indiv_spike_inds_s <- 1 + r + sort(sample(x = c(1:r.vec[s]), size = num_in_spike[2], replace = FALSE))
          }

          if (s != 1) {
            indiv_spike_inds_s <- 1 + r + sum(r.vec[1:(s-1)]) + sort(sample(x = c(1:r.vec[s]), size = num_in_spike[s], replace = FALSE))
          }
          spike_inds <- c(spike_inds, indiv_spike_inds_s)
        }

        gamma[spike_inds,] <- 0
      }

      diag(Sigma_beta)[gamma == 0] <- 1/1000
      beta <- matrix(list(), nrow = 1, ncol = 1)
      beta[[1,1]] <- matrix(MASS::mvrnorm(1, mu = rep(0, n_beta), Sigma = Sigma_beta), ncol = 1)

      # If desired, set spike coefficients to be identically 0
      if (identically_zero) {
        beta[[1,1]][gamma == 0] <- 0
      }
    }

    # Combine the Vs
    V.star.joint <- V
    if (r == 0) {
      V.star.joint[[1,1]] <- matrix(nrow = n, ncol = r)
    }

    Vs.star <- Vs
    for (s in 1:q) {
      if (r.vec[s] == 0) Vs.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
    }

    VStar <- cbind(1, do.call(cbind, V.star.joint), do.call(cbind, Vs.star))

    if (response == "binary") {
      Y <- EY <- matrix(list(), nrow = 1, ncol = 1)
      EY[[1,1]] <- pnorm(VStar %*% beta[[1,1]]) # True probability of being a case
      Y[[1,1]] <- matrix(rbinom(n, size = 1, prob = EY[[1,1]]), ncol = 1)
      tau2 <- matrix(list(), nrow = 1, ncol = 1)
      s2nY_coef <- NULL
    }

    if (response == "continuous") {

      # If a true error variance is provided, set tau2 to this value
      if (length(error_vars) > q) {
        tau2 <- matrix(list(), nrow = 1, ncol = 1)
        tau2[[1,1]] <- error_vars[q+1]
      }

      # If a true error variance is not provided, use prior
      if (length(error_vars) == q) {
        tau2 <- matrix(list(), nrow = 1, ncol = 1)
        tau2[[1,1]] <- matrix(1/rgamma(1, shape = shape, rate = rate))
      }

      Y <- EY <- matrix(list(), nrow = 1, ncol = 1)
      EY[[1,1]] <- VStar %*% beta[[1,1]]
      error_y <- matrix(rnorm(n, mean = 0, sd = sqrt(tau2[[1,1]])), ncol = 1)

      # -------------------------------------------------------------------------
      # Standardizing the variance of the signal in the response
      # -------------------------------------------------------------------------

      if (is.null(s2nY)) {
        s2nY_coef <- NULL
      }

      if (!is.null(s2nY)) {
        # Calculating the scaling coefficient so that the variance of the response = s2nY * noise variance
        s2nY_coef <- sqrt(s2nY) * sd(error_y)/sd(EY[[1,1]])
        EY[[1,1]] <- s2nY_coef * EY[[1,1]]
      }

      Y[[1,1]] <- EY[[1,1]] + error_y
    }
  }

  # -------------------------------------------------------------------------
  # Adding missingness if desired
  # -------------------------------------------------------------------------

  if (is.null(missingness)) {
    missing_data <- missing_obs <-  matrix(list(), nrow = q, ncol = 1)
    missing_obs_Y <- Y_missing <- matrix(list(), nrow = 1, ncol = 1)
  }

  if (!is.null(missingness)) {
    if (missingness != "missingness_in_data" & missingness != "both") { # If missing in response
      missing_data <- missing_obs <- matrix(list(), nrow = q, ncol = 1)
    }

    if (missingness != "missingness_in_response" & missingness != "both") { # If missing in data
      Y_missing <- missing_obs_Y <- matrix(list(), nrow = 1, ncol = 1)
    }

    if (missingness == "missingness_in_response" | missingness == "both") { # If missing in response
      missing_obs_Y <- matrix(list(), nrow = 1, ncol = 1)
      missing_obs_Y[[1,1]] <- sort(sample(1:n, size = prop_missing * n, replace = FALSE))

      Y_missing <- Y
      Y_missing[[1,1]][missing_obs_Y[[1,1]],] <- NA
    }

    if (missingness == "missingness_in_data" | missingness == "both") { # If missing in data

      missing_obs <- missing_data <- missing_cols <- missing_rows <- matrix(list(), nrow = q, ncol = 1)

      if (missing_data_type == "entrywise") { # if removing observations entrywise
        for (s in 1:q) {
          # these are counters going down the columns of R. So 9 would be the 9th entry counting down.
          missing_obs[[s,1]] <- sort(sample(x = 1:length(data[[s,1]]), size = prop_missing*length(data[[s,1]]), replace = FALSE))

          # Duplicate Xs so that I have one with the full data and one with the missing data
          missing_data[[s,1]] <- data[[s,1]]
          missing_data[[s,1]][missing_obs[[s,1]]] <- NA
        }
      }

      if (missing_data_type == "columnwise") { # if removing entire columns
        # Gives the column indices to remove

        for (s in 1:q) {
          # These are counters going down the COLUMNS of X. So 9 would be the 9th column.
          missing_cols[[s,1]] <- sort(sample(x=1:n, size = n*prop_missing, replace = FALSE))

          if (s != 1) {
            avail_obs <- c(1:n)[!(c(1:n) %in% unlist(missing_cols[1:(s-1),]))]
            missing_cols[[s,1]] <- sample(x=avail_obs, size = n*prop_missing, replace = FALSE)
          }

          # Duplicate Xs so that I have one with the full data and one with the missing data
          missing_data[[s,1]] <- data[[s,1]]
          missing_data[[s,1]][,missing_cols[[s,1]]] <- NA
          missing_obs[[s,1]] <- sort(which(is.na(missing_data[[s,1]])))
        }
      }

      if (missing_data_type == "MNAR") { # Restrict missing to within features
        for (s in 1:q) {
          # Sort the entries in source s
          sorted_entries_s <- sort(c(data[[s,1]]), decreasing = FALSE)

          # Calculate the threshold at which the bottom prop_missing are set to NA
          number_below_lod <- prop_missing * length(sorted_entries_s)
          LOD <- ceiling(max(sorted_entries_s[1:number_below_lod]))

          # Duplicate Xs so that I have one with the full data and one with the missing data
          missing_data[[s,1]] <- data[[s,1]]
          missing_data[[s,1]][missing_data[[s,1]] < LOD] <- NA
          missing_obs[[s,1]] <- sort(which(is.na(missing_data[[s,1]])))
        }
      }
    }
  }

  # -------------------------------------------------------------------------
  # Return
  # -------------------------------------------------------------------------

  list(data = data, # The "observed data"
       Y = Y, # The "observed outcome"
       missing_data = missing_data, # Missing data
       missing_obs = missing_obs, # Missing data
       Y_missing = Y_missing, missing_obs_Y = missing_obs_Y, # Missing data
       s2nX = s2nX, s2nX_coef = s2nX_coef, # Scaling for s2n in data
       s2nY = s2nY, s2nY_coef = s2nY_coef, # Scaling for the response
       joint.structure = joint.structure, # Joint structure
       indiv.structure = indiv.structure, # Individual structure
       overall.structure = overall.structure, # Joint + Individual structure
       V = V, U = U, Vs = Vs, W = W, # Components of the structure
       beta = beta, EY = EY, tau2 = tau2, gamma = gamma, p.prior = p.prior)
}


#' var_explained
#'
#' Calculates the proportion of variance explained in each data source by each
#' estimated structure.
#'
#' @param BSFP.fit Output from the BSFP function
#' @param iters_burnin (vector) If adding a burn-in, add the iterations after burn-in. May also include thinning.
#' Default is NULL if burn-in is already added.
#' @param source.names (vector) Vector of source names. Should be of length \eqn{q}, the number of sources.
#'
#' @details Calculate the proportion of variation in each data source by the estimated structures. Provides posterior
#' summaries of the variance explained.
#'
#' @returns Returns a list of posterior summaries for the variance explained by the joint and individual structures in each source.
#' \item{Joint}{Named list with each entry corresponding to the posterior mean and 95% credible interval for the proportion
#' of variance explained by the joint structure in each data source. }
#' \item{Individual}{Named list with each entry corresponding to the posterior mean and 95% credible interval for the proportion
#' of variance explained by the individual structures in the corresponding data source. }
#'
#' @export
#'
#' @examples
#' # Setting up the data
#' n <- 50
#' p.vec <- c(75, 100)
#' q <- 2
#'
#' # Setting up the model parameters
#' true_params <- list(error_vars = c(1,1),
#'                     joint_var = 1,
#'                     indiv_vars = c(1,1),
#'                     beta_vars = c(1, 1, rep(1, q)),
#'                     response_vars = c(shape = 1, rate = 1))
#'
#' # Choose ranks
#' r <- 3
#' r.vec <- c(3, 3)
#' ranks <- c(r, r.vec)
#'
#' # Number of posterior sampling iterations
#' nsample <- 1000
#' burnin <- nsample/2
#' iters_burnin <- (burnin+1):nsample
#'
#' # Generate data
#' data.c1 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE)
#'
#' # Run BSFP for 1000 iterations
#' bsfp.c1 <- bsfp(data = data.c1$data, Y = data.c1$Y, nsample = nsample)
#'
#' bsfp.c1.var.explained <- var_explained(BSFP.fit = bsfp.c1, iters_burnin = iters_burnin, source.names = c("Expression", "Methylation"))
var_explained <- function(BSFP.fit, iters_burnin = NULL, source.names) {

  # Save the standardized data
  data <- BSFP.fit$data

  # Save the number of sources
  q <- nrow(data)

  # If no iterations after burn-in given
  if (is.null(iters_burnin)) {
    iters_burnin <- 1:length(BSFP.fit$J.draw)
  }

  # Save the estimated joint and individual structures and the fitted Y values
  J.draw <- BSFP.fit$J.draw[iters_burnin]
  A.draw <- BSFP.fit$A.draw[iters_burnin]

  # Joint structure --
  joint_var_exp <- lapply(1:q, function(s) {
    var_by_iter <- sapply(1:length(iters_burnin), function(iter) {
      1 - (frob(data[[s,1]] - J.draw[[iter]][[s,1]]))/frob(data[[s,1]])
    })
  })
  names(joint_var_exp) <- source.names

  # Individual structure --
  indiv_var_exp <- lapply(1:q, function(s) {
    var_by_iter <- sapply(1:length(iters_burnin), function(iter) {
      1 - (frob(data[[s,1]] - A.draw[[iter]][[s,1]]))/frob(data[[s,1]])
    })
  })
  names(indiv_var_exp) <- source.names

  # Summarizing the results
  joint_summary <- lapply(1:q, function(s) {
    c(Mean = mean(joint_var_exp[[s]]), Lower = quantile(joint_var_exp[[s]], 0.025), Upper = quantile(joint_var_exp[[s]], 0.975))
  })
  indiv_summary <- lapply(1:q, function(s) {
    c(Mean = mean(indiv_var_exp[[s]]), Lower = quantile(indiv_var_exp[[s]], 0.025), Upper = quantile(indiv_var_exp[[s]], 0.975))
  })
  names(joint_summary) <- names(indiv_summary) <- source.names

  # Return
  list(Joint = joint_summary, Individual = indiv_summary)
}
