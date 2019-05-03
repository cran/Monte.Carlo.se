#' Standard Errors for Summaries Based on 2 or More Vectors of Monte Carlo Output
#'
#' \code{mc.se.matrix} --- jackknife and bootstrap SEs for statistics based on k correlated samples
#'
#' @param x N by k matrix or data frame of MC output, k>1
#' @param xcol vector specifying columns of x to use
#' @param B B=0 (default) means use jackknife, B>0 means use bootstrap with B resamples, If B>0, then a seed must be given to start the bootstrap resampling
#' @param seed seed=NULL (default) used with jackknife, otherwise needs positive integer
#' @param summary.f summary function with arguments (index,xdata), x goes into xdata
#' @param ... Additional arguments to be passed
#'
#' @return Returns data frame of summary.f, Monte Carlo (MC) SE of summary.f, MC sample size N,
#' method (jackknife or bootstrap), B and seed if bootstrap is used
#'
#' @details Suppose that k columns of Monte Carlo output are produced from N
#' independent Monte Carlo samples, and summary statistics involving 2 or more columns are to be
#' reported in tables.  \code{mc.se.matrix} gives Monte Carlo standard errors (SEs) for
#' these summary statistics (e.g., the ratio of 2 column means or variances). The vignette 
#' \code{vignette("Example2", package = "Monte.Carlo.se")} is a detailed 
#' account of using mc.se.vector.
#' 
#' @examples
#' \donttest{
#' # First create MC output used in Table 9.1, p. 367, of Boos and Stefanski (2013).
#' # norm15 holds 10,000 sample means, 20% trimmed means, and medians
#' # for normal samples of size 15
#'
#' N <- 10000
#' set.seed(346)                   # sets the random number seed
#' z <- matrix(rnorm(N*15),nrow=N) # N rows of N(0,1) samples, n=15
#' out.m.15 <- apply(z,1,mean)     # mean for each sample
#' out.t20.15 <- apply(z,1,mean,trim=0.20)   # 20% trimmed mean for each sample
#' out.med.15 <- apply(z,1,median)           # median for each sample
#'
#' # Save all 1000 blocks of 3 estimators in a data frame
#' norm15 <- data.frame(mean=out.m.15,trim20=out.t20.15,median=out.med.15)
#'
#' # Pearson correlation based on 2 vectors of MC output
#' # Note that the 2 vectors are in xdata, index is for the rows of xdata
#' corr<-function(index,xdata){cor(xdata[index,1],xdata[index,2])}
#'
#' # Compute jackknife SE for summary.f=corr
#'
#' mc.se.matrix(norm15,xcol=c(1,2),summary.f=corr)
#' #      summary          se     N    method
#' #  1 0.9367602 0.001256079 10000 Jackknife
#'
#' # Compute bootstrap SE for summary.f=corr
#'
#' mc.se.matrix(norm15,xcol=c(1,2),summary.f=corr,B=1000,seed=3928)
#' #      summary          se     N    method    B seed
#' #  1 0.9367602 0.001287065 10000 Bootstrap 1000 3928
#'
#' # Rerun with B=5000
#'
#' mc.se.matrix(norm15,xcol=c(1,2),summary.f=corr,B=5000,seed=3928)
#' #      summary          se     N    method    B seed
#' #  1 0.9367602 0.001266177 10000 Bootstrap 5000 3928
#'
#' # Compute jackknife SE for summary.f=ratio.var
#' # = ratio of variances of the two columns
#' # A ratio of 2 variances facilitates comparison of the variances
#'
#' ratio.var <- function(index,xdata)
#' {var(xdata[index,1])/var(xdata[index,2])}
#'
#' # ratio of column 1 variance to column 2 variance
#'
#' mc.se.matrix(norm15,xcol=c(1,2),summary.f=ratio.var)
#' #     summary          se     N    method
#' # 1 0.8895367 0.006263652 10000 Jackknife
#' # Coupled with SE=0.006, the ratio=0.89 shows the second variance is larger than the fitst
#'
#' # ratio of column 2 variance to column 1 variance
#' # Same conclusion as for the previous ratio
#'
#' mc.se.matrix(norm15,xcol=c(2,1),summary.f=ratio.var)
#' #    summary          se     N    method
#' # 1 1.124181 0.007915667 10000 Jackknife
#' }
#'
#' # Function Code
#'
#' mc.se.matrix <- function(x,xcol,B=0,seed=NULL,summary.f,...){
#'  # x is an n by k matrix or data frame of MC output, k>1
#'  # xcol is a vector specifying columns of x to use
#'  # summary.f is the summary function with arguments (index,xdata)
#'  # ... is for additional arguments to summary.f
#'  # B=0 means use jackkife, B>0 means use bootstrap with B resamples
#'  # If B>0, then a seed must be given to start the bootstrap resampling
#'  N=nrow(x)
#'  x=x[,xcol]   # use columns selected
#'  if(B>0){
#'    if(is.null(seed))stop('If B>0, then seed for the bootstrap must be given in the call')
#'    se=boot.se(1:N,B=B,theta = summary.f,xdata=x,...)}
#'  else{se=jack.se(1:N,theta = summary.f,xdata=x,...)}
#'  summ = summary.f(1:N,xdata=x,...)
#'  if(B>0) {out=data.frame(summary=summ,se,N,method="Bootstrap",B,seed)}
#'  else{out=data.frame(summary=summ,se,N,method="Jackknife")}
#'  return(out)}
#'  
#' @references
#' Boos, D. D., and Osborne, J. A. (2015), "Assessing Variability of Complex Descriptive
#' Statistics in Monte Carlo Studies using Resampling Methods," 
#' International Statistical Review, 25, 775-792.
#'
#' @export
#'
#' @author Dennis Boos, Kevin Matthew, Jason Osborne
#'
mc.se.matrix <- function(x,xcol,B=0,seed=NULL,summary.f,...){
  # x is an n by k matrix or data frame of MC output, k>1
  # xcol is a vector specifying columns of x to use
  # summary.f is the summary function with arguments (index,xdata)
  # ... is for additional arguments to summary.f
  # B=0 means use jackkife, B>0 means use bootstrap with B resamples
  # If B>0, then a seed must be given to start the bootstrap resampling
  N=nrow(x)
  x=x[,xcol]   # use columns selected
  if(B>0){
    if(is.null(seed))stop('If B>0, then seed for the bootstrap must be given in the call')
    se=boot.se(1:N,B=B,theta = summary.f,xdata=x,...)}
  else{se=jack.se(1:N,theta = summary.f,xdata=x,...)}
  summ = summary.f(1:N,xdata=x,...)
  if (B>0) {out=data.frame(summary=summ,se,N,method="Bootstrap",B,seed)}
  else{out=data.frame(summary=summ,se,N,method="Jackknife")}
  return(out)}
