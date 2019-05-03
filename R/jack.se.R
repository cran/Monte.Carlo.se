#' Jackknife Standard Error
#'
#' \code{jack.se} --- gives jackknife standard error (SE=estimated standard deviation)
#'
#' @param x vector of data
#' @param theta --- function (statistic) applied to the data (e.g., mean, median, var)
#' @param ... Additional arguments to be passed
#'
#' @return Returns the jackknife standard error of theta(x)
#'
#' @examples
#' \donttest{
#' # simple example, data from Boos and Osborne (2015, Table 3)
#' # using theta = coefficient of variation = mean/sd
#'
#' x=c(1,2,79,5,17,11,2,15,85)
#' cv=function(x){sd(x)/mean(x)}
#' cv(x)
#' # [1] 1.383577
#'
#' jack.se(x,theta=cv)
#' # [1] 0.3435321
#'
#' # More complex example using two samples, se for ratio of means
#' # data from Higgins (2003, problem 4.4, p. 142), LDH readings on 7 patients
#'
#' before=c(89,90,87,98,120,85,97)
#' after=c(76,101,84,86,105,84,93)
#'
#' # requires function using row index as "data,"
#' # real data is extra parameter xdata
#'
#' ratio.means <- function(index,xdata){
#'  mean(xdata[index,1])/mean(xdata[index,2])}
#'
#' # ratio of means for before-after data
#'
#' ratio.means(index=1:7,xdata=data.frame(before,after))
#' # [1] 1.058824
#'
#' # jackknife SE for ratio of means
#'
#' jack.se(x=1:7,theta=ratio.means,xdata=data.frame(before,after))
#' # [1] 0.03913484
#'
#' # To illustrate use with Monte Carlo output, first create some sample data
#' # 10,000 samples of size 15 from the Laplace (double exp) distribution
#'
#' N<-10000
#' set.seed(450)
#' z1 <- matrix(rexp(N*15),nrow=N)
#' z2 <- matrix(rexp(N*15),nrow=N)
#' z<-(z1-z2)/sqrt(2)              # subtract standard exponentials
#' out.m.15   <- apply(z,1,mean)
#' out.t20.15 <- apply(z,1,mean,trim=0.20)
#' out.med.15 <- apply(z,1,median)
#'
#' # The three datasets (out.m.15,out.t20.15,out.med.15) each contain 10,000 values.
#' # If we want use the variance of each column in a table, then to get
#' # the Monte Carlo standard error of those 3 variances,
#'
#' jack.se(out.m.15,theta = var)
#' # [1] 0.0009612314
#'
#' jack.se(out.t20.15,theta = var)
#' # [1] 0.0007008859
#'
#' jack.se(out.med.15,theta = var)
#' # [1] 0.0008130531
#' 
#' # Function Code
#' 
#' jack.se=function(x, theta, ...){
#'  call <- match.call()
#'  n <- length(x)
#'  u <- rep(0, n)
#'  for(i in 1:n) {u[i] <- theta(x[ - i], ...)}
#'  jack.se <- sqrt(((n - 1)/n) * sum((u - mean(u))^2))
#'  return(jack.se)}
#' }
#' @export
#'
#' @details The code was modified from code associated with the appendix of Efron and Tibshirani (1993)
#'
#' @references Efron and Tibshirani (1993), *An Introduction to the Bootstrap*.
#' 
#' Boos, D. D., and Osborne, J. A. (2015), "Assessing Variability of Complex Descriptive
#' Statistics in Monte Carlo Studies using Resampling Methods," 
#' International Statistical Review, 25, 775-792.
#'
#' @seealso \code{\link{boot.se}} --- \code{\link{mc.se.vector}} --- \code{\link{mc.se.matrix}}
#'
#' @author Dennis Boos, Kevin Matthew, Jason Osborne
#'
jack.se=function(x, theta, ...){
 call <- match.call()
 n <- length(x)
 u <- rep(0, n)
 for(i in 1:n) {u[i] <- theta(x[ - i], ...)}
 jack.se <- sqrt(((n - 1)/n) * sum((u - mean(u))^2))
 return(jack.se)}
