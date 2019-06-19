# Purpose: Package documentation
# Updated: 19/06/19

#' @useDynLib BNR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' BNR: Bivariate Normal Regression
#'
#' This package performs estimation and inference for the parameters of a
#' bivariate normal regression model. Estimation is performed using \code{\link{fit.bnr}}.
#' Inference on regression parameters for the target outcome is performed using
#' \code{\link{test.bnr}}.
#'
#' @author Zachary R. McCaw
#' @docType package
#' @name BNR
NULL
