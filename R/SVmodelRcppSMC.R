#' SVmodelRcppSMC: A package for example implementations of PMCMC inference
#'
#' This SVmodelRcppSMC package provides Particle Gibbs and Particle Gibbs with
#' ancestor sampling estimation of the standard SV model using the more general
#' RcppSMC package for the conditional SMC part of the algorithm.
#'
#'
#' 
## SVmodelRcppSMC namespace: start
#' @useDynLib SVmodelRcppSMC, .registration = TRUE
#' @importFrom stats rnorm
#' @importFrom rlang .data
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics legend matplot title
#' @importFrom stats dgamma dnorm rgamma runif
#' @importFrom MASS mvrnorm
## SVmodelRcppSMC namespace: end
#' @docType package
#' @name SVmodelRcppSMC
NULL
