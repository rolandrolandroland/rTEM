% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rTEM_features.R
\name{calc_T_val_observed}
\alias{calc_T_val_observed}
\title{T value for observed}
\usage{
calc_T_val_observed(
  observed,
  expected,
  r = NULL,
  func = "K",
  rmin = 0,
  rmax = NA,
  sqrt = FALSE,
  obs_name = c(NA, "cor_name"),
  expect_name = c("rrl_", "mmean"),
  K_cor = "trans",
  G_cor = "km",
  F_cor = "km",
  GXGH_cor = "km",
  GXHG_cor = "km"
)
}
\description{
T value for observed
}
