#' @importFrom foreach foreach %do% %dopar%
NULL

setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))

#' @export CopulaDistribution
CopulaDistribution <- setClass('CopulaDistribution',
                               slots = c(nu = 'numeric',
                                         gamma = 'numeric'),
                               prototype = list(nu = Inf))

#' @export CopulaDynamics
CopulaDynamics <- setClass('CopulaDynamics',
                           slots = c(alpha = 'numeric',
                                     beta = 'numeric',
                                     Omega = 'matrixOrNULL'),
                           prototype = list(alpha = 0, beta = 0, Omega = NULL))

#' @export CopulaSpecification
CopulaSpecification <- setClass('CopulaSpecification',
                                slots = c(distribution = 'CopulaDistribution',
                                          dynamics = 'CopulaDynamics'))

#' @export
copula_filter <- function(spec, u) {
  shocks <- .copula_shocks(spec, u)
  shocks_std <- .copula_shocks_std(spec, shocks)

  if (is.null(spec@dynamics@Omega)) {
    spec@dynamics@Omega <- .copula_Omega(spec, shocks_std)
  }

  Correlation <- .copula_Correlation(.copula_Q(spec, shocks_std))

  scores <- .copula_scores(spec, shocks, Correlation)

  list(
    spec = spec,
    shocks = shocks,
    shocks_std = shocks_std,
    Correlation = Correlation,
    scores = scores,
    ll = sum(scores)
  )
}

#' Simulate uniform residuals from the copula
#'
#' @param spec CopulaSpecification
#' @param n.sim Simulation horizon
#' @param m.sim Number of simulation
#' @param Q_initial What to set Q to initially (default: Omega)
#'
#' @return List of m.sim n.simxN uniform residuals
#' @export
copula_simulate <- function(spec, n.sim, m.sim, Q_initial = NULL) {
  if (is.null(Q_initial)) {
    stopifnot(!is.null(spec@dynamics@Omega))
    Q_initial <- spec@dynamics@Omega
  }

  uv_distributions <- .copula_uv_distributions(spec)

  foreach(i = 1:m.sim) %dopar% {
    shocks <- .copula_simulate_shocks(spec, n.sim, Q_initial)

    # Compute the uniform residuals using the marginal distributions of
    # the copula
    sapply(seq_along(uv_distributions),
           function(i) ghyp::pghyp(shocks[, i], uv_distributions[[i]]))
  }
}

#' @export
copula_is_constant <- function(spec) {
  all(c(spec@dynamics@alpha, spec@dynamics@beta) == 0)
}

.copula_rghyp <- function(n, spec, Correlation) {
  ghyp::rghyp(n, .copula_mv_distribution(spec, Correlation))
}


# Distributions ----------------------------------------------------------

.copula_uv_distributions <- function(spec) {
  lapply(spec@distribution@gamma, function(gamma_i) {
    if (is.infinite(spec@distribution@nu)) {
      return(ghyp::gauss())
    }

    ghyp::student.t(nu = spec@distribution@nu, gamma = gamma_i)
  })
}

.copula_mv_distribution <- function(spec, Correlation) {
  mu <- rep(0, ncol(Correlation))

  if (is.infinite(spec@distribution@nu)) {
    return(ghyp::gauss(mu = mu, sigma = Correlation))
  }

  ghyp::student.t(mu = mu, nu = spec@distribution@nu,
                  gamma = spec@distribution@gamma, sigma = Correlation)
}

# Filtering --------------------------------------------------------------

#' Compute shocks in the copula model
#'
#' @param spec A CopulaSpecification
#' @param u TxN matrix of uniforms
#'
#' @return shocks TxN matrix of distributed shocks
.copula_shocks <- function(spec, u) {
  uv_distributions <- .copula_uv_distributions(spec)

  shocks <- foreach(i = 1:ncol(u)) %dopar% {
    ghyp::qghyp(u[, i], uv_distributions[[i]], method = 'splines')
  }

  matrix(unlist(shocks), ncol = ncol(u))
}

#' Standardize shocks by forcing them to have expectation 0 and unit variance
#' and scale them by the conditional variance in the DCC model. See
#' Christoffersen for details and whatever paper he refers to for
#' justification.
#'
#' @param spec CopulaSpecification
#' @param shocks TxN shocks
#'
#' @return TxN matrix of standardized shocks
.copula_shocks_std <- function(spec, shocks) {
  uv_distributions <- .copula_uv_distributions(spec)
  T <- nrow(shocks)

  # Step 1: Standardize shocks in traditional sense of the word by subtracting
  # the distribution mean and dividing by standard deviation.
  standardize <- function(dist, i)
    (shocks[, i] - ghyp::mean(dist)) / sqrt(ghyp::vcov(dist))

  shocks <- rbind(mapply(standardize, uv_distributions, 1:ncol(shocks)))

  # Step 2: Following Christoffersen's first paper, we also scale the shocks
  # by the "conditional variance" for this period. Supposedly improves the
  # consistency of estimates.
  alpha <- spec@dynamics@alpha
  beta <- spec@dynamics@beta

  # Add an initial observation, assumed zero and set initial diagonal of Q = 1
  shocks <- rbind(0, shocks)
  shocks_std <- shocks * NA
  qdiag <- shocks * NA
  qdiag[1, ] <- 1

  for (t in 2:(T + 1)) {
    qdiag[t, ] <-
      (1 - alpha - beta) +
      beta * qdiag[t - 1, ] +
      alpha * (shocks[t - 1, ] / qdiag[t - 1, ]) ^ 2

    shocks_std[t, ] <- shocks[t, ] / sqrt(qdiag[t, ])
  }

  # Drop the initial observation
  shocks_std[-1, ]
}

.copula_Omega <- function(spec, shocks_std, use_cor = F) {
  # Christoffersen does write out the second formula, but this might be a
  # sloppy way of meaning the sample correlation; difference should be minimal
  # since shocks should be standardized, however, the sample might be "off"
  if (use_cor) {
    correlation <- cor(shocks_std)
  }
  else {
    correlation <- (t(shocks_std) %*% shocks_std) / nrow(shocks_std)
  }

  # TODO Add support for computing Gamma matrix here as well
  correlation
}

.copula_Q_tp1 <- function(spec, Q_t, shocks_std_t) {
  alpha <- spec@dynamics@alpha
  beta <- spec@dynamics@beta
  Omega <- spec@dynamics@Omega
  stopifnot(!is.null(Omega))

  (1 - alpha - beta) * Omega +
    beta * Q_t +
    alpha * shocks_std_t %*% t(shocks_std_t)
}

.copula_Q <- function(spec, shocks_std) {
  N <- ncol(shocks_std)
  T <- nrow(shocks_std)

  # One extra observation at the start
  Q <- array(dim = c(N, N, T + 1))
  shocks_std <- rbind(0, shocks_std)

  # t0 observation = unconditional; diag(1, N) is also plausible
  Q[,, 1] <- spec@dynamics@Omega

  for (t in 2:(T + 1)) {
    Q[,, t] <- .copula_Q_tp1(spec, Q[,, t - 1], shocks_std[t - 1, ])
  }

  Q[,, -1]
}

.copula_Correlation <- function(Q) {
  fn <- function(t) {
    Q_t <- Q[,, t]
    inv_sqrt <- solve(sqrt(diag(diag(Q_t))))
    inv_sqrt %*% Q_t %*% inv_sqrt
  }

  N <- dim(Q)[1]
  T <- dim(Q)[3]

  correlation_list <- lapply(seq(T), fn)

  array(unlist(correlation_list), dim = c(N, N, T))
}

# Simulation -------------------------------------------------------------

.copula_simulate_shocks <- function(spec, n.sim, Q_initial) {
  N <- ncol(Q_initial)

  shocks <- matrix(ncol = N, nrow = n.sim)
  shocks_std <- shocks
  Q <- array(dim = c(N, N, n.sim))
  Correlation <- Q

  # Initiate dynamics
  Q[,, 1] <- Q_initial
  Correlation[,, 1] <- .copula_Correlation(array(Q[,, 1], dim = c(N, N, 1)))

  # If the copula is constant, there is no need to update the dynamics.
  # We can just simulate n.sim shocks directly
  if (copula_is_constant(spec)) {
    shocks <- .copula_rghyp(n.sim, spec, Correlation[,, 1])
    return(shocks)
  }

  for (t in seq(n.sim)) {
    shocks[t, ] <- .copula_rghyp(1, spec, Correlation[,, t])

    if (t < n.sim) {
      shocks_std[t, ] <- .copula_shocks_std(spec, rbind(shocks[t, ]))
      Q[,, t + 1] <- .copula_Q_tp1(spec, Q[,, t], shocks_std[t, ])
      Correlation[,, t + 1] <- .copula_Correlation(array(Q[,, t + 1],
                                                         dim = c(N, N, 1)))
    }
  }

  shocks
}

# Log Likelihoods --------------------------------------------------------

.copula_scores <- function(spec, shocks, Correlation) {
  joint <- .copula_ll_joint(spec, shocks, Correlation)
  marginal <- .copula_ll_marginal(spec, shocks)

  joint - marginal
}

.copula_ll_marginal <- function(spec, shocks) {
  uv_distributions <- .copula_uv_distributions(spec)

  fn <- function(dist, i) ghyp::dghyp(shocks[, i], dist, logvalue = TRUE)

  ll <- foreach(i = 1:ncol(shocks), .combine = 'cbind') %do% {
    ghyp::dghyp(shocks[, i], uv_distributions[[i]], logvalue = TRUE)
  }

  cbind(rowSums(ll))
}

.copula_ll_joint <- function(spec, shocks, Correlation) {
  fn <- function(t) {
    mvtnorm::dmvnorm(shocks[t, ], sigma = Correlation[,, t], log = TRUE)
  }

  # Use our special density function that cuts of a lot of the fat
  if (is.finite(spec@distribution@nu)) {
    fn <- function(t) {
      temp <- .dghst(rbind(shocks[t, ]),
              nu = spec@distribution@nu,
              gamma = spec@distribution@gamma,
              sigma = Correlation[,, t])

      unname(as.vector(temp))
    }
  }

  # If you ever feel the need to verify our own density function, uncomment
  # the below lines to use ghyp package. It takes CONSIDERABLY longer.
  # fn <- function(t) {
  #   mv_distribution <- .copula_mv_distribution(spec, Correlation[,, t])
  #   ghyp::dghyp(shocks[t, ], mv_distribution, logvalue = TRUE)
  # }

  # Running this in parallel actually seems to have little benefit
  ll <- lapply(seq(nrow(shocks)), fn)
  cbind(unlist(ll))
}
