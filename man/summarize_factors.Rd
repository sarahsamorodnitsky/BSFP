% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bsfp.R
\name{summarize_factors}
\alias{summarize_factors}
\title{summarize_factors}
\usage{
summarize_factors(
  data,
  Y = NULL,
  iters_burnin,
  aligned_results,
  ranks,
  tau2.draw = NULL,
  Xm.draw = NULL,
  Ym.draw = NULL
)
}
\arguments{
\item{data}{A matrix of lists or a list of matrices that share the same number of
columns. The matrices must be oriented in pxn orientation. May contain NAs if
there are missing values in the dataset.}

\item{Y}{A matrix of lists or a nx1 matrix of continuous or binary outcome.
May be NULL if no outcome is given. May contain NAs if there are missing outcomes.}

\item{iters_burnin}{(vector): indices for posterior samples after burn-in}

\item{aligned_results}{(list): results from match_align_bsfp}

\item{ranks}{Estimated joint and individual ranks from BSFP}

\item{tau2.draw}{Posterior samples for error variance in Y if response is given}

\item{Xm.draw}{Imputed values for missing values in X}

\item{Ym.draw}{Imputed values for unobserved outcomes in Y}
}
\description{
Calculate posterior summaries of aligned estimated factors from BSFP.
}