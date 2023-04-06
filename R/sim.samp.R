#' Code for Generating Random Samples
#'
#' \code{sim.samp} --- generates N samples of size n from distribution DIST
#'
#' @param nrep = number of random data sets
#' @param n = sample size of each data set
#' @param DIST = distribution to generate from, e.g., runif, rnorm, additional parameters in ...
#' @param ... Additional arguments to be passed
#'
#' @return N by n matrix, each row is a data set of size n
#'
#' @author Dennis Boos, Kevin Matthew, Jason Osborne 
#'
#' @importFrom stats var sd cor
#'
#' @examples
#' N<-100
#'  set.seed(346)            # sets the random number seed
#'  sim.samp(N,15,rnorm)->z  # 100 N(0,1) samples of size n=15
#'  sim.samp(N,40,rnorm,mean=10,sd=5)->z  # 100 N(10,25) samples of size n=40
#'
#' # Function Code
#'
#' sim.samp <- function(nrep,n,DIST,...){
#`  # simulates nrep samples from DIST of size n
#'  data <- matrix(DIST(n * nrep, ...), ncol = n, nrow = nrep)}
#'
#' @export
#'
sim.samp <- function(nrep,n,DIST,...){
  # simulates nrep samples from DIST of size n
  data <- matrix(DIST(n * nrep, ...), ncol = n, nrow = nrep)}


