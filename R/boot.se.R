#' Bootstrap Standard Error
#'
#' \code{boot.se} -- gives a bootstrap standard error (SE=estimated standard deviation)
#'
#' @param x vector of data
#' @param B number of bootstrap resamples (replications)
#' @param theta function (statistic) applied to the data (e.g., mean, median, var)
#' @param ... Additional arguments to be passed
#'
#' @return Returns the bootstrap SE of theta(x)
#'
#' @importFrom stats var sd cor
#'
#' @examples
#' \donttest{
#' # simple example, data from Boos and Osborne (2105, Table 3)
#' # using theta=coefficient of variation mean/sd
#'
#' x=c(1,2,79,5,17,11,2,15,85)
#' cv=function(x){sd(x)/mean(x)}
#'
#' cv(x)
#' # [1] 1.383577
#'
#' # bootstrap SE using B=1000
#'
#' set.seed(384)
#' boot.se(x,B=1000,theta=cv)
#' # [1] 0.3416897
#'
#' # More complex example using two samples, se for ratio of means
#' # data from Higgins (2003, problem 4.4, p. 142), LDH readings on 7 patients
#'
#' before=c(89,90,87,98,120,85,97)
#' after=c(76,101,84,86,105,84,93)
#'
#' # requires function using row index as "data", real data is extra parameter xdata
#' ratio.means <- function(index,xdata)
#'  {mean(xdata[index,1])/mean(xdata[index,2])}
#'
#' ratio.means(index=1:7,xdata=data.frame(before,after))
#' # [1] 1.058824
#'
#' # boostrap SE for ratio of means
#'
#' set.seed(2917)
#' boot.se(x=1:7,B=1000,theta=ratio.means,xdata=data.frame(before,after))
#' # [1] 0.03576659
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
#'
#' # The three datasets (out.m.15,out.t20.15,out.med.15) each contain 10000 values.
#' # If we want use the variance of each column in a table, then to get
#' # the Monte Carlo standard error of those 3 variances,
#'
#' set.seed(250)
#' boot.se(out.m.15,B=1000,theta = var)
#' # [1] 0.0009373835
#'
#' boot.se(out.t20.15,B=1000,theta = var)
#' # [1] 0.0007086057
#'
#' boot.se(out.med.15,B=1000,theta = var)
#' # [1] 0.0008307258
#' } # ends donttest
#' # Function Code
#' 
#' boot.se<-function(x, B, theta, ...){
#'   call <- match.call()
#'   n <- length(x)
#'   bootsam <- matrix(sample(x, size = n * B, replace = T), nrow = B)
#'   thetastar <- apply(bootsam, 1, theta, ...)
#'   se <- sd(thetastar)
#'   return(se)
#'   }
#'
#' @export
#'
#' @references Efron and Tibshirani (1993), _An Introduction to the Bootstrap_.
#' 
#' Boos, D. D., and Osborne, J. A. (2015), "Assessing Variability of Complex Descriptive
#' Statistics in Monte Carlo Studies using Resampling Methods," 
#' International Statistical Review, 25, 775-792.
#'
#' @seealso \code{\link{jack.se}} -- \code{\link{mc.se.vector}} --
#'  \code{\link{mc.se.matrix}}
#'
#' @details The code was modified from code associated with the appendix of Efron and Tibshirani (1993)
#'
#' @author Dennis Boos, Kevin Matthew, Jason Osborne
#'
boot.se<-function(x, B, theta, ...){
  call <- match.call()
  n <- length(x)
  bootsam <- matrix(sample(x, size = n * B, replace = T), nrow = B)
  thetastar <- apply(bootsam, 1, theta, ...)
  se <- sd(thetastar)
  return(se)}
