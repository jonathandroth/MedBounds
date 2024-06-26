% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_sharp_null_fsst_binary_m.R
\name{test_sharp_null_fsst_binary_m}
\alias{test_sharp_null_fsst_binary_m}
\title{Hypothesis test for the sharp null}
\usage{
test_sharp_null_fsst_binary_m(
  df,
  d,
  m,
  y,
  ordering = NULL,
  B = 500,
  cluster = NULL,
  weight.matrix = "diag",
  ats_only = F,
  alpha = 0.05,
  kappa = alpha/10,
  use_hybrid = T,
  num_Ybins = NULL,
  analytic_variance = FALSE,
  lambda = "dd"
)
}
\arguments{
\item{df}{A data frame}

\item{d}{Name of the treatment variable in the df}

\item{m}{Name of the mediator variable}

\item{y}{Name of the outcome variable, which is assumed to take a discrete
support}

\item{ordering}{A list with length equal to the cardinality of the support of
the mediator variable. The name of each element corresponds to a point in
the support, and each element is a vector that collects all m values that
are less than or equal to this point. If ordering = NULL, the standard
ordering is used.}

\item{B}{Bootstrap size, default is zero}

\item{cluster}{Cluster for bootstrap}

\item{weight.matrix}{Weight matrix used to implement FSST. Possible options
are "diag", "avar", "identity." Defaults is "diag" as in FSST.}

\item{num_Ybins}{(Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)}
}
\description{
This function tests the sharp null of Y(1,m) = Y(0,m). The
  outcome and mediator are both assumed to take finitely many different
  values. The inference is via applying Andrews, Roth, and Pakes
}
