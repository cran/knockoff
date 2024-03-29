% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats_sqrt_lasso.R
\name{stat.sqrt_lasso}
\alias{stat.sqrt_lasso}
\title{Importance statistics based on the square-root lasso}
\usage{
stat.sqrt_lasso(X, X_k, y, ...)
}
\arguments{
\item{X}{n-by-p matrix of original variables.}

\item{X_k}{n-by-p matrix of knockoff variables.}

\item{y}{vector of length n, containing the response variables of numeric type.}

\item{...}{additional arguments specific to \code{slim}.}
}
\value{
A vector of statistics \eqn{W} of length p.
}
\description{
Computes the signed maximum statistic
  \deqn{W_j = \max(Z_j, \tilde{Z}_j) \cdot \mathrm{sgn}(Z_j - \tilde{Z}_j),}
where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of 
\eqn{\lambda} at which the jth variable and its knockoff, respectively,
enter the SQRT lasso model.
}
\details{
With default parameters, this function uses the package \code{RPtests}
to run the SQRT lasso. By specifying the appropriate optional parameters, 
one can use different Lasso variants including Dantzig Selector, LAD Lasso,
SQRT Lasso and Lq Lasso for estimating high dimensional sparse linear models.

For a complete list of the available additional arguments, see \code{\link[RPtests]{sqrt_lasso}}.
}
\examples{
set.seed(2022)
p=50; n=50; k=10
mu = rep(0,p); Sigma = diag(p)
X = matrix(rnorm(n*p),n)
nonzero = sample(p, k)
beta = 3.5 * (1:p \%in\% nonzero)
y = X \%*\% beta + rnorm(n)
knockoffs = function(X) create.gaussian(X, mu, Sigma)

# Basic usage with default arguments
result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=stat.sqrt_lasso)
print(result$selected)

# Advanced usage with custom arguments
foo = stat.sqrt_lasso
k_stat = function(X, X_k, y) foo(X, X_k, y, q=0.5)
result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
print(result$selected)

}
\seealso{
Other statistics: 
\code{\link{stat.forward_selection}()},
\code{\link{stat.glmnet_coefdiff}()},
\code{\link{stat.glmnet_lambdadiff}()},
\code{\link{stat.lasso_coefdiff_bin}()},
\code{\link{stat.lasso_coefdiff}()},
\code{\link{stat.lasso_lambdadiff_bin}()},
\code{\link{stat.lasso_lambdadiff}()},
\code{\link{stat.random_forest}()},
\code{\link{stat.stability_selection}()}
}
\concept{statistics}
