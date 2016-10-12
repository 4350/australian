#' @importFrom foreach foreach
NULL

setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))

CopulaDistribution <- setClass('CopulaDistribution',
                               slots = c(nu = 'numeric',
                                         gamma = 'numeric'),
                               prototype = list(nu = Inf))

CopulaDynamics <- setClass('CopulaDynamics',
                           slots = c(alpha = 'numeric',
                                     beta = 'numeric',
                                     Omega = 'matrixOrNULL'),
                           prototype = list(alpha = 0, beta = 0, Omega = NULL))

CopulaSpecification <- setClass('CopulaSpecification',
                                slots = c(distribution = 'CopulaDistribution',
                                          dynamics = 'CopulaDynamics'))

copula_filter <- function(spec, u) {
  tic('shocks')
  shocks <- .copula_shocks(spec, u)
  toc()
  tic('shocks_std')
  shocks_std <- .copula_shocks_std(spec, shocks)
  toc()

  tic('Omega')
  if (is.null(spec@dynamics@Omega)) {
    spec@dynamics@Omega <- .copula_Omega(spec, shocks_std)
  }
  toc()

  tic("Correlation")
  Correlation <- .copula_Correlation(.copula_Q(spec, shocks_std))
  toc()
  tic("Scores")
  scores <- .copula_scores(spec, shocks, Correlation)
  toc()

  list(
    spec = spec,
    shocks = shocks,
    shocks_std = shocks_std,
    Correlation = Correlation,
    scores = scores,
    ll = sum(scores)
  )
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
#' @param spec
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
#' @param spec
#' @param shocks
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

.copula_Q <- function(spec, shocks_std) {
  N <- ncol(shocks_std)
  T <- nrow(shocks_std)
  alpha <- spec@dynamics@alpha
  beta <- spec@dynamics@beta
  Omega <- spec@dynamics@Omega

  stopifnot(!is.null(Omega))

  # One extra observation at the start
  Q <- array(dim = c(N, N, T + 1))
  shocks_std <- rbind(0, shocks_std)

  # t0 observation = unconditional; diag(1, N) is also plausible
  Q[,, 1] <- Omega

  for (t in 2:(T + 1)) {
    Q[,, t] <-
      (1 - alpha - beta) * Omega +
      beta * Q[,, t - 1] +
      alpha * (shocks_std[t - 1, ] %*% t(shocks_std[t - 1, ]))
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

# Log Likelihoods --------------------------------------------------------

.copula_scores <- function(spec, shocks, Correlation) {
  tic('joint')
  joint <- .copula_ll_joint(spec, shocks, Correlation)
  toc()
  tic('marginal')
  marginal <- .copula_ll_marginal(spec, shocks)
  toc()

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
    ghyp::dghyp(shocks[t, ], .copula_mv_distribution(spec, Correlation[,, t]),
                logvalue = TRUE)
  }

  if (!all(spec@distribution@gamma == 0)) {
    stopifnot(is.finite(spec@distribution@nu))

    fn <- function(t) {
      temp <- .dghst(rbind(shocks[t, ]),
              nu = spec@distribution@nu,
              gamma = spec@distribution@gamma,
              sigma = Correlation[,, t])

      unname(as.vector(temp))
    }
  }

  cbind(unlist(
    foreach(t = 1:nrow(shocks),
            .export = c('.copula_mv_distribution', '.dghst')) %dopar% fn(t)))
}
