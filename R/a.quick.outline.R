#' Monte.Carlo.se: A package for computing standard errors of Monte Carlo estimates
#'
#' The Monte.Carlo.se package has two main R functions
#' \code{\link{mc.se.vector}} and \code{\link{mc.se.matrix}} that are built on two basic functions
#' \code{\link{jack.se}} and \code{\link{boot.se}} for computing jackknife and bootstrap standard errors, respectively.
#'
#' A Monte Carlo Study often results in an N by k matrix X of estimates based on N independent samples of
#' simulation data.  X might contain parameter estimates (like the mean or median) or test results
#' (say 0 for "accept the null hypothesis" or 1 for "reject the null"). Summary calculations from X,
#' such as means or variances, are then displayed in a table or plot.
#'
#' The purpose of the Monte.Carlo.se package is to compute estimates of the standard deviation
#' (called standard errors) of these table entries.
#'
#' Examples are included in the documentation for \code{\link{mc.se.vector}} and \code{\link{mc.se.matrix}}.
#' However, more extended examples may be found in the three vignettes found by clicking "index" at
#' the bottom of this page and then on "User guides, package vignettes and other documentation."
#' 
#' @author Dennis Boos, Kevin Matthew, Jason Osborne
#'
#' @name A.Quick.Outline
NULL

