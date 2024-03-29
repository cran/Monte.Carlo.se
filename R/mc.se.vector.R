#' Standard Errors for Monte Carlo Output Summaries
#'
#' \code{mc.se.vector} --- gives jackknife and bootstrap SEs for a vector of MC output (say N estimates)
#'
#' @param x vector from N independent Monte Carlo replications
#' @param summary.f summary function computed from x (e.g., mean, median, var)
#' @param B B=0 means use jackknife (default), B>0 means use bootstrap with B resamples, If B>0, then a seed must be given to start the bootstrap resampling
#' @param seed seed=NULL (default) used with jackknife, otherwise needs a positive integer
#' @param ... Additional arguments to be passed
#'
#' @return Returns data frame of summary.f , MC standard error of summary.f, MC sample size N,
#' method (jackknife or bootstrap), B and seed if bootstrap is used
#'
#' @details Suppose that an N-vector of Monte Carlo output (thus, a sample of size N)
#' is produced from N
#' independent Monte Carlo samples, and a summary statistic like the mean or variance is to be
#' reported in a table.
#' \code{mc.se.vector} gives Monte Carlo standard errors (SEs) for these summary statistics.
#' The vignette \code{vignette("Example1", package = "Monte.Carlo.se")} is a detailed account of using \code{mc.se.vector}.
#'
#' @importFrom stats var sd cor
#'
#' @examples
#' \donttest{
#' # First create MC output used for Table 9.1, p. 367, of Boos and Stefanski (2013)
#' # norm15 holds 10,000 sample means, 20% trimmed means, and medians
#' # computed from normal samples of size 15
#'
#' N<-10000
#' set.seed(346)                   # sets the random number seed
#' z <- matrix(rnorm(N*15),nrow=N) # N rows of N(0,1) samples, n=15
#' out.m.15 <- apply(z,1,mean)     # mean for each sample
#' out.t20.15 <- apply(z,1,mean,trim=0.20)   # 20% trimmed mean for each sample
#' out.med.15 <- apply(z,1,median)           # median for each sample
#'
#' # Save all 1000 blocks of 3 estimators in a data frame
#'
#' norm15 <- data.frame(mean=out.m.15,trim20=out.t20.15,median=out.med.15)
#'
#' # Compute the jackknife Standard Error (SE) for summary.f=mean
#' # This summary is useful for estimating the bias of estimators.
#'
#' mc.se.vector(norm15[,1],summary.f=mean)
#' #        summary         se     N    method
#' # 1 -0.001920642 0.00256714 10000 Jackknife
#'
#' # Compute a bootstrap SE for summary.f=mean
#'
#' mc.se.vector(norm15[,1],B=1000,seed=4822,summary.f=mean)
#' #        summary          se     n    method    B seed
#' # 1 -0.001920642 0.002573516 10000 Bootstrap 1000 4822
#'
#' # compare to basic R
#'
#' mean(norm15[,1])
#' # [1] -0.001920642
#'
#' sd(norm15[,1])/sqrt(10000)
#' # [1] 0.00256714
#'
#' # Illustrate use with summaries having additional parameters.
#' # Compute the jackknife SE for summary.f=varn
#' # Multiplying the variance by sample size n allows comparison for different n
#'
#' varn=function(x,n){n*var(x)}
#'
#' # The additional parameter n replaces the ... in the mc.se.vector definition
#'
#' mc.se.vector(norm15[,1],summary.f=varn,n=15)
#' #     summary         se     N    method
#' # 1 0.9885314 0.01407249 10000 Jackknife
#'
#' # Compute a bootstrap SE for summary.f=varn
#'
#' mc.se.vector(norm15[,1],B=1000,seed=3029,summary.f=varn,n=15)
#' #     summary         se     n    method    B seed
#' # 1 0.9885314 0.01408173 10000 Bootstrap 1000 3029
#' }
#'
#' # Function Code
#' 
#' mc.se.vector <- function(x,summary.f,B=0,seed=NULL,...){
#' # x is a vector from N independent Monte Carlo replications
#' # summary.f is a summary function computed from MC output
#' # ... is for additional arguments to summary.f
#' # B=0 means use jackkife, B>0 means use bootstrap with B resamples
#'  if (B>0){
#'   if(is.null(seed))stop('If B>0, then seed for the bootstrap must be given in the call')
#'   set.seed(seed);se=boot.se(x,B=B,theta = summary.f,...)}
#'  else{se = jack.se(x,theta = summary.f,...)}
#' summ = summary.f(x,...)
#' if(B>0){out=data.frame(summary=summ,se,n=length(x),method="Bootstrap",B,seed)}
#'  else{out=data.frame(summary=summ,se,N=length(x),method="Jackknife")}
#' return(out)}
#' 
#' @references 
#' Boos, D. D., and Osborne, J. A. (2015), "Assessing Variability of Complex Descriptive
#' Statistics in Monte Carlo Studies using Resampling Methods," 
#' International Statistical Review, 25, 775-792.
#'
#' @author Dennis Boos, Kevin Matthew, Jason Osborne
#'
#' @export
#'
mc.se.vector <- function(x,summary.f,B=0,seed=NULL,...){
  # x is a vector from N independent Monte Carlo replications
  # summary.f is a summary function computed from MC output
  # ... is for additional arguments to summary.f
  # B=0 means use jackkife, B>0 means use bootstrap with B resamples
  # If B>0, then a seed must be given to start the bootstrap resampling
  if (B>0) {
    if(is.null(seed))stop('If B>0, then seed for the bootstrap must be given in the call')
    set.seed(seed);se=boot.se(x,B=B,theta = summary.f,...)}
  else{se = jack.se(x,theta = summary.f,...)}
  summ = summary.f(x,...)
  if (B>0) {out=data.frame(summary=summ,se,n=length(x),method="Bootstrap",B,seed)}
  else{out=data.frame(summary=summ,se,N=length(x),method="Jackknife")}
  return(out)}

