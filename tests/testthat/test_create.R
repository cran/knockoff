test_that('Equicorrelated knockoffs have the right correlation structure', {
  n = 20; p = 10
  X = normc(rnorm_matrix(n,p))
  X_ko_default = knockoff.create(X, 'equicorrelated', randomized=F)
  X_ko_randomized = knockoff.create(X, 'equicorrelated', randomized=T)
  
  G = t(X) %*% X
  s = min(2*min(eigen(G)$values), 1)
  for (X_ko in list(X_ko_default, X_ko_randomized)) {
    expect_equal(t(X_ko) %*% X_ko, G)
    expect_equal(t(X) %*% X_ko, G - diag(s,p,p))
  }
})

test_that('SDP knockoffs have the right correlation structure', {
  skip_on_cran()
  if (!has_cvxpy())
    skip('CVXPY not available')
  
  n = 20; p = 10
  X = normc(rnorm_matrix(n,p))
  X_ko_default = knockoff.create(X, 'sdp', randomized=F)
  X_ko_randomized = knockoff.create(X, 'sdp', randomized=T)
  
  offdiag <- function(A) A - diag(diag(A))
  G = t(X) %*% X
  tol = 1e-4
  for (X_ko in list(X_ko_default, X_ko_randomized)) {
    expect_equal(t(X_ko) %*% X_ko, G, tolerance=tol)
    expect_equal(offdiag(t(X) %*% X_ko), offdiag(G), tolerance=tol)
    expect_true(all(diag(t(X) %*% X_ko) < 1+tol))
  }
})