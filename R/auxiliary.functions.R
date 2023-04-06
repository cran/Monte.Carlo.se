#' Auxiliary Functions
#'
#' @description Auxiliary functions to be used in the Monte.Carlo.se package,
#' mainly with mc.se.matrix.  Scroll down to the Examples Section to
#' see the actual code.
#'
#' @param index index = usually of the form 1:N
#' @param xdata actual data
#' @param x Input vector in calls to jack.var and boot.var
#' @param theta theta = function in calls to jack.var and boot.var
#' @param B Bootstrap reps in calls to boot.var
#' @param ... Additional arguments to be passed
#' @param true true parameter value when computing mean squared error
#' @param n sample size
#'
#' @importFrom stats var sd cor
#'
#' @examples
#' # These are extra functions included in the MCse package
#' # The following functions are to be used with mc.se.matrix
#'
#' ratio.var <- function(index,xdata)         # ratio of variances
#'  {var(xdata[index,1])/var(xdata[index,2])}
#'
#' # The above function is for the ratio of the sample variance of column 1 to
#' # the sample variance  of column 2 of xdata.
#' # Note that the actual data goes into xdata, the second argument of ratio.var.
#' # Example call for 10,000 means and medians:
#' # ratio.var(1:10000,xdata=cbind(out.m.15,out.med.15))
#'
#' ratio.sd<-function(index,xdata){           # ratio of standard deviations
#'  sd(xdata[index,1])/sd(xdata[index,2])}
#'
#' ratio.mse<-function(index,xdata,true){     # ratio of mean squared errors
#'  mean((xdata[index,1]-true)^2)/mean((xdata[index,2]-true)^2)}
#'
#' ratio.mean.vhat.var<-function(index,xdata){# estimates in col 1, vhats in col. 2
#'  mean(xdata[index,2])/var(xdata[index,1])}
#'
#' ratio.mean.sdhat.sd<-function(index,xdata){# estimates in col 1, SEs in col. 2
#'  mean(xdata[index,2])/sd(xdata[index,1])}
#'
#' corr<-function(index,xdata){               # simple correlation
#'  cor(xdata[index,1],xdata[index,2])}
#'
#' # These next two functions correspond to jack.se and boot.se.
#' # x is a data vector, and theta is a function applied to x.
#' # Each returns a variance estimate for theta(x).
#'
#' jack.var <- function(x, theta, ...){  # jackknife estimate of variance
#'  n <- length(x)
#'  u <- rep(0, n)
#'  for(i in 1:n){u[i] <- theta(x[ - i], ...)}
#'  jack.var <-((n-1)/n)* sum((u-mean(u))^2)
#'  return(jack.var)}
#'
#' boot.var <- function(x,B,theta, ...){ # bootstrap estimate of variance
#'  n <- length(x)
#'  bootsam <- matrix(sample(x,size = n*B,replace=T), nrow=B)
#'  thetastar <- apply(bootsam,1,theta,...)
#'  boot.var <- var(thetastar)
#'  return(boot.var)}
#'
#' @author Dennis Boos, Kevin Matthew, Jason Osborne
#' @name Extra
NULL
#' @rdname Extra
ratio.var <- function(index,xdata)         # ratio of variances
  {var(xdata[index,1])/var(xdata[index,2])}
#' @rdname Extra
ratio.sd<-function(index,xdata){           # ratio of standard deviations
   sd(xdata[index,1])/sd(xdata[index,2])}
#' @rdname Extra
ratio.mse<-function(index,xdata,true){    # ratio of mean squared errors
  mean((xdata[index,1]-true)^2)/mean((xdata[index,2]-true)^2)}
#' @rdname Extra
ratio.mean.vhat.var<-function(index,xdata){   # estimates in col 1, vhats in col. 2
  mean(xdata[index,2])/var(xdata[index,1])}
#' @rdname Extra
ratio.mean.sdhat.sd<-function(index,xdata){   # estimates in col 1, SEs in col. 2
  mean(xdata[index,2])/sd(xdata[index,1])}
#' @rdname Extra
corr<-function(index,xdata){                  # simple correlation
  cor(xdata[index,1],xdata[index,2])}
#' @rdname Extra
cv <- function(x){sd(x)/mean(x)}              # coefficient of variation
#' @rdname Extra
varn<-function(x,n){n*var(x)}                 # sample variance times sample size n
#' @rdname Extra
jack.var <- function(x, theta, ...){          # jackknife estimate of variance
  n <- length(x)
  u <- rep(0, n)
  for(i in 1:n){u[i] <- theta(x[ - i], ...)}
  jack.var <-((n-1)/n)* sum((u-mean(u))^2)
  return(jack.var)}
#' @rdname Extra
boot.var <- function(x,B,theta, ...){         # bootstrap estimate of variance
  n <- length(x)
  bootsam <- matrix(sample(x,size = n*B,replace=T), nrow=B)
  thetastar <- apply(bootsam,1,theta,...)
  boot.var <- var(thetastar)
  return(boot.var)}
