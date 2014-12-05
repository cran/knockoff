#' Create knockoff variables
#' 
#' Creates knockoff variables for the original variables.
#' 
#' @param X normalized n-by-p design matrix (n >= 2p)
#' @param method either 'equicorrelated' or 'sdp'
#' @param randomized whether the knockoffs are deterministic or randomized
#' @return The n-by-p knockoff matrix
#' 
#' @export
knockoff.create <- function(X, method=c('equicorrelated','sdp'), randomized=F) {
  fn = switch(match.arg(method), 
              equicorrelated = create_equicorrelated,
              sdp = create_sdp)
  fn(X, randomized)
}

# Create equicorrelated knockoffs.
create_equicorrelated <- function(X, randomized) {
  # Compute SVD and U_perp.
  X.svd = decompose(X, randomized)
  
  # Set s = min(2 * smallest eigenvalue of X'X, 1), so that all the correlations
  # have the same value 1-s.
  if (any(X.svd$d <= 1e-5 * max(X.svd$d)))
    stop(paste('Data matrix is rank deficient.',
               'Equicorrelated knockoffs will have no power.'))
  lambda_min = min(X.svd$d)^2
  s = min(2*lambda_min, 1)
  
  # Construct the knockoff according to Equation 1.4.
  X_ko = (X.svd$u %*diag% (X.svd$d - s / X.svd$d) +
          X.svd$u_perp %*diag% (sqrt(2*s - (s/X.svd$d)^2))) %*% t(X.svd$v)
}

# Create SDP knockoffs.
create_sdp <- function(X, randomized) {
  if (!has_cvxpy())
    stop('To use SDP knockoffs, you must have Python and the CVXPY library.')
  
  # Compute SVD and U_perp.
  X.svd = decompose(X, randomized)
  
  # Check for rank deficiency.
  tol = 1e-5
  d = X.svd$d
  d_inv = 1 / d
  d_zeros = d <= tol*max(d)
  if (any(d_zeros)) {
    warning(paste('Data matrix is rank deficient.',
                  'Model is not identifiable, but proceeding with SDP knockoffs'))
    d_inv[d_zeros] = 0
  }
  
  # Compute the Gram matrix and its (pseudo)inverse.
  G = (X.svd$v %*diag% d^2) %*% t(X.svd$v)
  G_inv = (X.svd$v %*diag% d_inv^2) %*% t(X.svd$v)
  
  # Optimize the parameter s of Equation 1.3 using SDP.
  s = solve_sdp(G)
  s[s <= tol] = 0
  
  # Construct the knockoff according to Equation 1.4:
  C.svd = canonical_svd(2*diag(s) - (s %diag*% G_inv %*diag% s))
  X_ko = X - (X %*% G_inv %*diag% s) + 
    (X.svd$u_perp %*diag% sqrt(pmax(0, C.svd$d))) %*% t(C.svd$v)
}

# Compute the SVD of X and construct an orthogonal matrix U_perp such that
# U_perp * U = 0.
decompose <- function(X, randomized) {
  n = nrow(X); p = ncol(X)
  stopifnot(n >= 2*p)
  
  result = canonical_svd(X)
  Q = qr.Q(qr(cbind(result$u, matrix(0,n,p))))
  u_perp = Q[,(p+1):(2*p)]
  if (randomized) {
      Q = qr.Q(qr(rnorm_matrix(p,p)))
      u_perp = u_perp %*% Q
  }
  result$u_perp = u_perp
  result
}

# Solves the semidefinite programming problem:
#
#   maximize    sum(s)
#   subject to  0 <= s <= 1
#               2G - diag(s) >= 0
#
# Because R lacks a decent SDP library, we call out to Python. It would be nice
# to use an R-to-Python FFI, but there is only one maintained package (rPython)
# and it suffers from several defects:
#
#   1) No MS Windows support
#
#   2) The Python interpreter is fixed at compile time. Consequently any OS X
#      user not using the system Python (i.e., every OS X user) will have to
#      compile rPython from scratch (which is not as easy as one might expect).
#
# So we call the Python interpreter directly, using JSON as the data 
# serialization format.
solve_sdp <- function(G) {
  source_file = system.file('python', 'solve_sdp.py', package='knockoff')
  G.json = toJSON(G, collapse=' ')
  s.json = system2('python', source_file, stdout=T, input=G.json)
  if (!is.null(attr(s.json, 'status')))
    stop('Error calling Python SDP solver.')
  s = fromJSON(s.json)
  return(s)
}