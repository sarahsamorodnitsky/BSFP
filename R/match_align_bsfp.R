#' match_align_bsfp
#'
#' Correcting for rotation, permutation, and sign invariance by applying the
#' MatchAlign algorithm (Poworoznek et al. (2021)) to results from BSFP fit.
#'
#' @param BSFP.fit (list): a list of results from the BPMF model fit
#' @param y (matrix of lists): outcome vector if available
#' @param model_params (list): model parameters used in model fitting
#' @param p.vec (vector): vector with number of variables in each source
#' @param iters_burnin (vector): indices for posterior samples after burn-in
#' @param piv.list (list of matrices or ints): list of pivot matrices for joint then
#' individual structures or indices for pivots from posterior samples
#' @param index (int): for sensitivity of the results, choose a pivot that falls
#' the some indices away from the chosen pivot (index specifies the number of indices)
#'
#' @details The estimated loadings and scores from BSFP are not identifiable due to
#' rotation, permutation, and sign invariance. This invariance obstructs using posterior
#' summaries to study the estimated factors. To address these issues, we adapt the MatchAlign
#' algorithm proposed by Poworoznek et al. (2021) to address these ambiguities. This function
#' formats the posterior samples for the loadings and scores in the appropriate format to
#' run the MatchAlign alignment and returns the final aligned factors. The factors within each structure
#' are ordered from most-to-least variance explained.
#'
#' This also allows users to attempt
#' different pivots than the one used by default to examine sensitivity.
#'
#' @references
#' \insertAllCited{}
#'
#' @returns This function returns the aligned loadings, scores, and regression coefficients.
#' \item{joint.scores.final}{List of joint scores after alignment.}
#' \item{joint.loadings.final}{List of joint loadings after alignment.}
#' \item{joint.betas.final}{List of regression coefficients corresponding to joint factors
#' after alignment.}
#' \item{joint_pivot_index}{Returns the index after burn-in of the pivot used in the alignment.
#' This can be used to further examine sensitivity.}
#' \item{individual.scores.final}{List of individual scores for each source after alignment.}
#' \item{individual.loadings.final}{List of individual loadings for each source after alignment.}
#' \item{individual.betas.final}{List of regression coefficients corresponding to individual
#' factors after alignment.}
#' \item{individual_pivot_index}{Returns the indices after burn-in of the pivot used in the alignment
#' of the individual factors for each source. This can be used to further examine sensitivity.}
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

match_align_bsfp <- function(BSFP.fit, y = NULL, model_params, p.vec, iters_burnin, piv.list = NULL, index = 0) {

  # Save the ranks from the model fit
  ranks <- BSFP.fit$ranks
  joint.rank <- ranks[1]
  indiv.ranks <- ranks[-1]

  # Save the number of datasets
  q <- nrow(BSFP.fit$data)

  # Save the number of posterior samples
  nsample <- length(BSFP.fit$V.draw)

  # Save the burnin
  burnin <- length(iters_burnin)

  # If no response is given
  if (is.null(y)) {
    joint.betas.final <- individual.betas.final <- joint_var_betas <- NULL
    indiv_var_betas <- lapply(1:q, function(s) NULL)
  }

  # Create a vector for the indices of each rank
  if (!is.null(y)) {

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

    beta.ind <- lapply(1:(q+1), function(s) {

      if (is.null(rank.inds[[s]])) {
        NULL
      }

      else {
        rank.inds[[s]] + 1
      }
    })
  }

  # Checking if a pivot was provided
  if (is.null(piv.list)) {
    piv.list <- lapply(1:(q+1), function(i) NULL)
  }

  # Joint structure (plus joint regression coefficients) --

  # Save the overall number of features
  p <- sum(p.vec)

  # Combine joint loadings and betas (fix betas at 0 if rank is 0)
  if (!is.null(y)) {

    if (r > 0) {
      joint.loadings <- lapply(iters_burnin, function(iter) {
        rbind(do.call(rbind, BSFP.fit$U.draw[[iter]]), # Joint loadings
              t(BSFP.fit$beta.draw[[iter]][[1,1]][beta.ind[[1]],])) # Joint regression coefficients
      })
    }

    if (r == 0) {
      joint.loadings <- lapply(iters_burnin, function(iter) {
        rbind(do.call(rbind, BSFP.fit$U.draw[[iter]]), # Joint loadings
              t(matrix(0, nrow = 1, ncol = 1))) # Joint regression coefficients
      })
    }

  }

  if (is.null(y)) {
    joint.loadings <- lapply(iters_burnin, function(iter) {
      rbind(do.call(rbind, BSFP.fit$U.draw[[iter]]))
    })
  }

  # Save joint scores
  joint.scores <- lapply(iters_burnin, function(iter) BSFP.fit$V.draw[[iter]][[1,1]])

  # Save the prior variance on the betas for joint factors
  if (!is.null(y)) {
    joint_var_betas <- diag(rep(model_params$beta_vars[2], ranks[1]))
  }

  # If we estimated more than 1 factor
  if (joint.rank > 1) {
    # Apply the factor switching method to the joint structure
    joint.results.rotate <- jointRot_multi(joint.loadings, joint.scores, y = y, var_betas = joint_var_betas, index = index, piv = piv.list[[1]])

    # Separate the joint loadings from the joint scores
    joint.loadings.final <- lapply(joint.results.rotate$lambda, function(iter) iter[1:p,])

    if (!is.null(y)) {
      joint.betas.final <- lapply(joint.results.rotate$lambda, function(iter) t(iter[p+1,,drop=FALSE]))
    }

    # Applying the permutation and sign switching to the scores
    joint.scores.final <- joint.results.rotate$eta

    # Save the pivot index
    joint_pivot_index <- joint.results.rotate$pivot_index
  }

  # If we estimated 1 or fewer factors
  if (joint.rank <= 1) {
    joint.scores.final <- joint.scores
    joint.loadings.final <- lapply(joint.loadings, function(iter) iter[1:p,,drop=FALSE])

    if (!is.null(y)) {
      joint.betas.final <- lapply(joint.loadings, function(iter) t(iter[p+1,,drop=FALSE]))
    }

    joint_pivot_index <- NULL
  }

  # Individual structure (plus individual regression coefficients) --

  # Combine the individual loadings and betas
  if (!is.null(y)) {
    individual.loadings <- lapply(1:q, function(s) {

      if (indiv.ranks[s] > 0) {
        lapply(iters_burnin, function(iter) {
          rbind(BSFP.fit$W.draw[[iter]][[s,s]],
                BSFP.fit$beta.draw[[iter]][[1,1]][beta.ind[[s+1]],])
        })
      } else {
        lapply(iters_burnin, function(iter) {
          rbind(BSFP.fit$W.draw[[iter]][[s,s]],
                matrix(0, nrow = 1, ncol = 1))
        })
      }

    })
  }

  if (is.null(y)) {
    individual.loadings <- lapply(1:q, function(s) lapply(iters_burnin, function(iter) {
      rbind(BSFP.fit$W.draw[[iter]][[s,s]])
    }))
  }

  individual.scores <- lapply(1:q, function(s) lapply(iters_burnin, function(iter) BSFP.fit$Vs.draw[[iter]][[1,s]]))

  # Save the prior variance on the betas for the individual factors
  if (!is.null(y)) {
    indiv_var_betas <- lapply(1:q, function(s) diag(rep(model_params$beta_vars[s+2], ranks[s+1])))
  }

  # Applying the alignment algorithm
  individual.results.rotate <- lapply(1:q, function(s) {

    if (indiv.ranks[s] > 1) {
      out <- jointRot_multi(individual.loadings[[s]], individual.scores[[s]], y = y, var_betas = indiv_var_betas[[s]], index = index, piv = piv.list[[s+1]])
      out
    } else {
      list(lambda = individual.loadings[[s]], eta = individual.scores[[s]])
    }
  })

  # Save the final individual loadings and betas
  individual.loadings.final <- lapply(1:q, function(s) lapply(individual.results.rotate[[s]]$lambda, function(iter) iter[1:p.vec[s],,drop=FALSE]))

  if (!is.null(y)) {
    individual.betas.final <- lapply(1:q, function(s) lapply(individual.results.rotate[[s]]$lambda, function(iter) t(iter[p.vec[s]+1,,drop=FALSE])))
  }

  # Applying the permutation and sign switching to the scores
  individual.scores.final <- lapply(1:q, function(s) individual.results.rotate[[s]]$eta)

  # Save the pivot index
  individual_pivot_index <- lapply(1:q, function(s) {

    if (indiv.ranks[s] > 1) {
      individual.results.rotate[[s]]$pivot_index
    }

  })

  # ---------------------------------------------------------------------------
  # Reordering aligned results so that factors are ordered from most-to-least
  # variance explained.
  # ---------------------------------------------------------------------------

  # Order the components by squared Frobenius-norm of the corresponding rank-1 structure after burn-in --

  # For each component, calculate the posterior mean of the corresponding rank-1 structure
  if (!is.null(y)) {

    if (joint.rank > 0) {
      joint.rank1.structure <- lapply(1:joint.rank, function(r) {
        # Calculate the rank-1 structure for the given component at each iteration
        factor_r_structure <- lapply(1:burnin, function(iter) {
          rbind(joint.loadings.final[[iter]][,r,drop=FALSE], t(joint.betas.final[[iter]][r,,drop=FALSE])) %*%
            t(joint.scores.final[[iter]][,r,drop=FALSE])
        })

        # Calculate the posterior mean
        Reduce("+", factor_r_structure)/length(factor_r_structure)
      })
    }

  }

  if (is.null(y)) {

    if (joint.rank > 0) {
      joint.rank1.structure <- lapply(1:joint.rank, function(r) {
        # Calculate the rank-1 structure for the given component at each iteration
        factor_r_structure <- lapply(1:burnin, function(iter) {
          joint.loadings.final[[iter]][,r,drop=FALSE] %*%
            t(joint.scores.final[[iter]][,r,drop=FALSE])
        })

        # Calculate the posterior mean
        Reduce("+", factor_r_structure)/length(factor_r_structure)
      })
    }

  }

  # Calculate the norm of each rank-1 structure
  if (joint.rank > 0) {
    joint.structure.norm <- sapply(joint.rank1.structure, function(str) frob(str))
  }

  if (joint.rank == 0) {
    joint.structure.norm <- 1
  }

  # Order the factors
  joint.factor.order <- order(joint.structure.norm, decreasing = TRUE)

  # Reorder the joint scores, loadings, and betas after burn-in
  joint.scores.final.order <- lapply(1:burnin, function(iter) {
    joint.scores.final[[iter]][,joint.factor.order,drop=FALSE]
  })
  joint.loadings.final.order <- lapply(1:burnin, function(iter) {
    joint.loadings.final[[iter]][,joint.factor.order,drop=FALSE]
  })

  if (!is.null(y)) {
    joint.betas.final.order <- lapply(1:burnin, function(iter) {
      joint.betas.final[[iter]][joint.factor.order,,drop=FALSE]
    })
  }

  # Order the factors in each individual structure by the Frobenius norm of their corresponding rank-1 structure --

  # For each source and each component, calculate the posterior mean of the corresponding rank-1 structure
  if (!is.null(y)) {

    individual.rank1.structure <- lapply(1:q, function(s) {

      if (indiv.ranks[s] > 0) {
        lapply(1:indiv.ranks[s], function(r) {
          # Calculate the rank-1 structure for the given component at each iteration
          factor_r_structure <- lapply(1:burnin, function(iter) {
            rbind(individual.loadings.final[[s]][[iter]][,r,drop=FALSE], t(individual.betas.final[[s]][[iter]][r,,drop=FALSE])) %*%
              t(individual.scores.final[[s]][[iter]][,r,drop=FALSE])
          })

          # Calculate the posterior mean
          Reduce("+", factor_r_structure)/length(factor_r_structure)
        })
      } else {
        NULL
      }

    })
  }

  if (is.null(y)) {

    individual.rank1.structure <- lapply(1:q, function(s) {

      if (indiv.ranks[s] > 0) {
        lapply(1:indiv.ranks[s], function(r) {
          # Calculate the rank-1 structure for the given component at each iteration
          factor_r_structure <- lapply(1:burnin, function(iter) {
            individual.loadings.final[[s]][[iter]][,r,drop=FALSE] %*%
              t(individual.scores.final[[s]][[iter]][,r,drop=FALSE])
          })

          # Calculate the posterior mean
          Reduce("+", factor_r_structure)/length(factor_r_structure)
        })
      } else {
        NULL
      }

    })
  }

  # Calculate the norm of each rank-1 structure for each source
  individual.structure.norm <- lapply(1:q, function(s) {
    if (indiv.ranks[s] > 0) {
      sapply(individual.rank1.structure[[s]], function(str) frob(str))
    } else {
      1
    }
  })

  # Order the factors
  individual.factor.order <- lapply(1:q, function(s) {
    order(individual.structure.norm[[s]], decreasing = TRUE)
  })

  # Reorder the individual scores, loadings, and betas after burn-in
  individual.scores.final.order <- lapply(1:q, function(s) {
    lapply(1:burnin, function(iter) {
      individual.scores.final[[s]][[iter]][,individual.factor.order[[s]],drop=FALSE]
    })
  })
  individual.loadings.final.order <- lapply(1:q, function(s) {
    lapply(1:burnin, function(iter) {
      individual.loadings.final[[s]][[iter]][,individual.factor.order[[s]],drop=FALSE]
    })
  })

  if (!is.null(y)) {
    individual.betas.final.order <- lapply(1:q, function(s) {
      lapply(1:burnin, function(iter) {
        individual.betas.final[[s]][[iter]][individual.factor.order[[s]],,drop=FALSE]
      })
    })
  }

  # If no response was given
  if (is.null(y)) {
    joint.betas.final.order <- individual.final.burnin <- NULL
  }

  # Return the final loadings, scores, and betas
  list(joint.scores.final = joint.scores.final.order,
       joint.loadings.final = joint.loadings.final.order,
       joint.betas.final = joint.betas.final.order,
       joint_pivot_index = joint_pivot_index,
       individual.scores.final = individual.scores.final.order,
       individual.loadings.final = individual.loadings.final.order,
       individual.betas.final = individual.betas.final.order,
       individual_pivot_index = individual_pivot_index)
}

#' MatchAlign Algorithm (code adapted from infinitefactor package by Evan Poworoznek (2021))
#'
#' Given results from BSFP, apply Varimax rotation and greedy alignment
#' using the MatchAlign algorithm to adjust for rotational, permutation,
#' and sign invariance. This function is used within the match_align_bsfp
#' function only. Code was
#' adapted from the \code{infinitefactor} package developed by Evan Poworoznek.
#'
#' @param lambda (list): list of loadings matrices. For BSFP, this will be combined
#'                loadings and betas.
#' @param eta (list): list of scores matrices
#' @param piv (matrix OR int): user-provided pivot. If NULL, use median largest singular value.
#'                      If matrix, use as pivot. If integer, use corresponding Varimax-rotated posterior sample
#' @param index (int): allows user to choose a different pivot for sensitivity analysis.
#'                     0 uses the default pivot. Could be positive or negative.
#'
#' @details This is an interval function that runs the MatchAlign algorithm,
#' originally developed by Evan Poworoznek, on structures estimated by the BSFP model.
#' Allows users to specify their own pivot or use a pivot several indices away from the default pivot.
#' As is done in Poworoznek et al. (2021), this function uses the median condition number to identify a pivot.
#'
#' @references
#' \insertAllCited{}
#'
#' @returns This function returns a list of aligned loadings and scores, as well as the index of the
#' pivot used in the alignment.
#' \item{lambda}{List of loadings after alignment.}
#' \item{eta}{List of scores after alignment.}
#' \item{pivot_index}{Int specifying the index of the pivot used in the alignment.}
#'
#' @export

jointRot_multi <- function(lambda, eta, piv = NULL, index = 0) {

  # Apply the Varimax rotation to the loadings
  vari = lapply(lambda, varimax)

  # Save the resulting Varimax-rotated loadings
  loads = lapply(vari, `[[`, 1)

  # Save the resulting rotation matrix
  rots = lapply(vari, `[[`, 2)

  # Apply the rotation matrix to the corresponding scores
  rotfact = mapply(`%*%`, eta, rots, SIMPLIFY = FALSE)

  # If no user-provided pivot is given, use default
  if (is.null(piv)) {
    norms = sapply(loads, norm, "2")
    piv = loads[order(norms)][[round(length(lambda)/2) + index]]

    # Save the posterior sample after burn-in that was used as pivot
    piv_index_to_return <- c(1:length(loads))[order(norms)][round(length(lambda)/2) + index]
  }

  # If user-provided pivot is given
  if (is.matrix(piv)) {
    piv_index_to_return <- NULL
  }

  if (!is.matrix(piv)) {
    piv = loads[[piv]]

    # Save the posterior sample after burn-in that was used as pivot
    piv_index_to_return <- piv
  }

  # Match to the defined pivot
  matches = lapply(loads, infinitefactor::msfOUT, piv)
  lamout = mapply(infinitefactor::aplr, loads, matches, SIMPLIFY = FALSE)
  etaout = mapply(infinitefactor::aplr, rotfact, matches, SIMPLIFY = FALSE)

  # Return
  return(list(lambda = lamout, eta = etaout, pivot_index = piv_index_to_return))
}
