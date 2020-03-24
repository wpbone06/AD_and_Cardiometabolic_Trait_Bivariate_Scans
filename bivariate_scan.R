## Date/Time : 07/22/2016
## Author    : Wei Zhao (@Einstein @UPenn)
## Support   : Benjamin Voight (@UPenn) and Danish Saleheen (@UPenn)
## R function: df2Test (Version 1.0)
##    This function is used to perform the 2-degree-freedom test under a bivariate normal distribution.
## Parameters:
##    x: a vector of z scores with length of two
##    mu: the estimated vector of means of the bivariate normal distribution
##    sigma.inv: the inverse of the covariance matrix of the bivariate normal distribution
##               It can be obtained by applying the solve() function on the covariance matrix.
## Returns   :
##    The -log10(pvalue) of this test.  Some times the p-values retruned can be quite small,
##    and R would approximate it to zero.  To avoid such situations, the -log10(pvalue) will
##    be derived and returned.
## Usage     :
##    To incorporate the function into an R script, one can use the source() function.
##    For example, if this df2Test.R is located at "/home/document/df2Test.R", then in
##    an R script, one can use source("/home/document/df2Test.R") to bring in the function.
##    This function can be used with the apply() function to perform the 2-degree-freedom test on many
##    vectors with z-scores.

df2Test <- function(x, mu, sigma.inv) {

        z <- x - mu
        stat <- t(z) %*% sigma.inv %*% z
        mlp <- - pchisq(q=stat, df=2, log.p=TRUE, lower.tail=FALSE) * log10(exp(1))
        return(mlp)

}
