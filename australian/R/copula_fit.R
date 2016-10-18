
#' Get theta and constraints for a copula model
#'
#' @param N how many factors?
#' @param distribution norm, t or ghst
#' @param constant constant copula?
#' @param upsilon with upsilon?
#'
#' @return list with theta, ui and ci
#' @export
copula_fit_theta <- function(N, distribution = 'norm', constant = T, upsilon = F) {
  # Distribution boundaries
  INITIAL_NU <- 8
  UI_NU <- rbind(1, -1)
  CI_NU <- rbind(6, -20)
  GAMMA_B <- 0.25

  # Initial and restrictions on alpha and beta
  INITIAL_ALPHABETA <- c(0.06, 0.91)
  UI_ALPHABETA <- rbind(diag(1, 2), c(-1, -1))
  CI_ALPHABETA <- rbind(0, 0, -0.9999)

  # Initial value and restrictions on phi (with Upsilon)
  INITIAL_PHI <- 0
  UI_PHI <- rbind(1, -1)
  CI_PHI <- rbind(0, -1)

  # Only support a single theta
  INITIAL_THETA <- 0

  pars <- list()
  ui <- list()
  ci <- list()

  # Distribution parameters to optimize over
  if (distribution == 't') {
    pars$distribution <- c(INITIAL_NU)
    names(pars$distribution) <- c('nu')

    ui$distribution <- UI_NU
    ci$distribution <- CI_NU
  }
  else if (distribution == 'ghst') {
    pars$distribution <- c(INITIAL_NU, rep(0, N))
    names(pars$distribution) <- c('nu', paste0('gamma', 1:N))

    # Add gamma restrictions
    ui$distribution <- magic::adiag(UI_NU,
                                    rbind(diag(1, N), diag(-1, N)))
    ci$distribution <- rbind(CI_NU,
                             cbind(rep(-GAMMA_B, N)), # lower
                             cbind(rep(-GAMMA_B, N))) # uppper
  }

  if (!constant) {
    pars$alphabeta <- INITIAL_ALPHABETA
    names(pars$alphabeta) <- c('alpha', 'beta')

    ui$alphabeta <- UI_ALPHABETA
    ci$alphabeta <- CI_ALPHABETA
  }

  if (upsilon) {
    pars$upsilon <- c(INITIAL_PHI, INITIAL_THETA)
    names(pars$upsilon) <- c('phi', 'theta')

    # theta is unrestricted
    ui$upsilon <- cbind(UI_PHI, 0)
    ci$upsilon <- CI_PHI
  }

  if (length(ui)) {
    ui <- do.call(magic::adiag, ui)
  }
  else {
    ui <- NULL
  }

  list(
    pars = unlist(pars),
    ui = ui,
    ci = unlist(ci)
  )
}
