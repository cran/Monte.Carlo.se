#' Standard Errors for Paired Comparisons of Monte Carlo Output Summaries
#'
#' \code{pairwise.se} --- gives jackknife and bootstrap SEs for all the pairwise difference or
#' ratios of Monte Carlo summaries.
#'
#' @param x vector from N independent Monte Carlo replications
#' @param xcol columns of x to be used
#' @param diff If TRUE (default), uses differences; if diff=F, uses ratios
#' @param digits number of digits to retain in output data frame
#' @param summary.f summary function computed from x (e.g., mean, median, var)
#' @param B B=0 means use jackknife (default), B>0 means use bootstrap with B resamples, 
#' If B>0, then a seed must be given to start the bootstrap resampling
#' @param seed seed=NULL (default) used with jackknife, otherwise needs a positive integer
#' @param ... Additional arguments to be passed
#'
#' @return Returns a data frame of the indiviual ith and jth column summaries (summi and summj),
#' the differences or ratios of those summaries (summary), MC standard error of the 
#' difference or ratio, MC sample size N,
#' method (jackknife or bootstrap), B and seed if bootstrap is used
#'
#' @details Suppose that an N-vector of Monte Carlo output (thus, a sample of size N)
#' is produced from N
#' independent Monte Carlo samples, and a summary statistic like the mean or variance is to be
#' reported in a table.
#' \code{pairwise.se} gives Monte Carlo standard errors (SEs) for all pairwise differences
#' or ratios of these summary statistics.
#' The vignette \code{vignette("Example3", package = "Monte.Carlo.se")} 
#' is a detailed account of using \code{pairwise.se}.
#'
#' @importFrom stats var sd cor
#'
#' @examples
#' \donttest{
#' # Using the output data matrix hold generated in vignette Example3, 
#' # calculate jackknife and bootstrap standard errors
#' # for the differences and ratios of the CV estimates.
#'
#' # Jackknife SE of Differences of CVs
#' 
#' # pairwise.se(hold,xcol=10:12,summary.f=cv)
#' # elem  summi  summj summary     se  t.stat    N    method
#' # 1 10 11 0.6884 0.7030 -0.0146 0.0299 -0.4877 1000 Jackknife
#' # 2 10 12 0.6884 0.6489  0.0395 0.0195  2.0274 1000 Jackknife
#' # 3 11 12 0.7030 0.6489  0.0541 0.0311  1.7374 1000 Jackknife
#' 
#' # Jackknife SE of Ratios of CVs
#' 
#' # pairwise.se(hold,xcol=10:12,diff=FALSE,summary.f=cv)
#' # elem  summi  summj summary     se  t.stat    N    method
#' # 1 10 11 0.6884 0.7030  0.9792 0.0429 -0.4833 1000 Jackknife
#' # 2 10 12 0.6884 0.6489  1.0608 0.0321  1.8972 1000 Jackknife
#' # 3 11 12 0.7030 0.6489  1.0833 0.0475  1.7531 1000 Jackknife
#' 
#' # Bootstrap SE of Differences of CVs
#' 
#' # pairwise.se(hold,xcol=10:12,B=1000,seed=770,summary.f=cv)
#' # elem  summi  summj summary     se  t.stat    B seed    N    method
#' # 1 10 11 0.6884 0.7030 -0.0146 0.0278 -0.5250 1000  770 1000 Bootstrap
#' # 2 10 12 0.6884 0.6489  0.0395 0.0182  2.1671 1000  770 1000 Bootstrap
#' # 3 11 12 0.7030 0.6489  0.0541 0.0303  1.7844 1000  770 1000 Bootstrap
#' 
#' # Bootstrap SE of Ratios of CVs
#' 
#' # pairwise.se(hold,xcol=10:12,diff=FALSE,B=1000,seed=770,summary.f=cv)
#' # elem  summi  summj summary     se  t.stat    B seed    N    method
#' # 1 10 11 0.6884 0.7030  0.9792 0.0390 -0.5316 1000  770 1000 Bootstrap
#' # 2 10 12 0.6884 0.6489  1.0608 0.0292  2.0797 1000  770 1000 Bootstrap
#' # 3 11 12 0.7030 0.6489  1.0833 0.0430  1.9372 1000  770 1000 Bootstrap
#' }
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
pairwise.se <- function(x,xcol,diff=TRUE,digits=4,B=0,seed=NULL,summary.f,...){
  # Returns the standard error of all pairwise differences or ratios
  # of summaries of x computed from columns of x specified by xcol
  # if diff=TRUE, differences are used. If diff=FALSE, then ratios.
  # x is an N by k matrix or data frame of MC output, k>1
  # summary.f is a simple summary function like mean, var, cv
  # ... is for additional arguments to summary.f
  # B=0 means use jackkife, B>0 means use bootstrap with B resamples
  # If B>0, then a seed must be given to start the bootstrap resampling
  # digits is the number of digits to retain
N=nrow(x)
k=length(xcol)
if(diff){mu=0} else{mu=1}
if(diff){combo=function(index,xdata,...){
  summary.f(xdata[index,1],...)-summary.f(xdata[index,2],...)}}
else{combo=function(index,xdata,...){
  summary.f(xdata[index,1],...)/summary.f(xdata[index,2],...)}}
if (B>0){if(is.null(seed))stop('If B>0, then seed for the bootstrap must be given in the call')}
for(i in 1:(k-1)){
  for(j in (i+1):k){
    if(B>0){out.df=mc.se.matrix(x=x[,c(xcol[i],xcol[j])],xcol=c(1,2),B=B,seed=seed,summary.f=combo,...)}
    else{out.df=mc.se.matrix(x=x[,c(xcol[i],xcol[j])],xcol=c(1,2),summary.f=combo,...)}
    summi = summary.f(x[,xcol[i]],...)
    summj = summary.f(x[,xcol[j]],...)
    elem=paste(xcol[i],xcol[j])
    out.temp=data.frame(elem,summi,summj,out.df)
    if(i==1 & j==2){out=out.temp}
    else{out=rbind(out,out.temp)}
  } # ends i loop
}  # ends j loop
if(B==0){out=data.frame(elem=out$elem,summi=round(out$summi,digits),
                        summj=round(out$summj,digits),summary=round(out$summary,digits),
                        se=round(out$se,digits),t.stat=round((out$summary-mu)/out$se,digits),N=out$N,method=out$method)}
else{out=data.frame(elem=out$elem,summi=round(out$summi,digits),
                    summj=round(out$summj,digits),summary=round(out$summary,digits),
                    se=round(out$se,digits),t.stat=round((out$summary-mu)/out$se,digits),
                    B=out$B,seed=out$seed,N=out$N,method=out$method)}
return(out)}


