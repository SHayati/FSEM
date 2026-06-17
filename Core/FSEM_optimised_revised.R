## =====================================================================
## em.estimation.optimised.revised
## ---------------------------------------------------------------------
## Wrapper around `em.estimation.optimised` that fixes an EM degeneracy
## for SEM models with coupled latents (e.g. main_simulation_3.R).
##
## Bug being fixed:
##   `initial.param` initialises every gamma.param entry to exactly 0.
##   With gamma = 0 the structural matrix B is zero, so the SEM-implied
##   prior covariance of the stacked latents Sigma_eta is block diagonal.
##   Since each indicator loads on a single latent, the loading matrix A
##   is also block diagonal, so the joint posterior
##       p(eta | z) = N( Sigma_eta A' Sigma_z^{-1} z , ... )
##   factorises across latents.  The eta samples drawn in the E-step are
##   therefore posterior-independent, the M-step OLS of eta1 on eta2
##   returns ~0, and EM is stuck at gamma = 0 forever.
##   Empirically: with true |gamma(t)| ~ 1.3, the original optimised
##   estimator returns |gamma(t)| ~ 0.09 (>>10x bias) and the other
##   parameters (lambda, phi, nu, ...) inherit the bias because the
##   sampled eta are mis-correlated.
##
## Fix:
##   Inject a small i.i.d. Gaussian jitter into every gamma.param entry
##   *before* calling em.estimation.optimised.  This breaks the
##   block-diagonal symmetry so the first M-step picks up a real signal
##   and gamma can grow.  Identification/masking is unaffected because
##   the structural diagMatrix factors zero out inactive entries inside
##   the algorithm.
##
## Extra arguments (everything else is forwarded verbatim):
##   gamma.init.sd : standard deviation of the gamma jitter.
##                   Default 0.3 (about a quarter of the typical gamma
##                   scale used in main_simulation_3.R).  Set to 0 to
##                   recover the original behaviour.
##   seed          : optional integer; if non-NULL, used to seed the
##                   gamma perturbation for reproducibility.  The EM's
##                   internal RNG state is restored afterwards.
##
## Important: this file MUST be sourced AFTER Core/FSEM_optimised.R.
## It does not modify Core/FSEM.R or Core/FSEM_optimised.R in place.
## =====================================================================

## =====================================================================
## .fsem_effective_df
## ---------------------------------------------------------------------
## Effective degrees of freedom (EDF) and an EDF-based modified AIC for a
## fitted FSEM, computed entirely from the quantities the EM stores in the
## returned object (`fit$hist`) plus the final variance parameters.
##
## Rationale.  Both the factor loadings and the structural coefficients are
## estimated in the M-step by a penalized (ridge-type) least squares of the
## form
##        coef = (A + delta * P)^{-1} b ,
## where A is the (Monte-Carlo averaged) Gram matrix t(F) Sigma^{-1} F, P is
## the roughness penalty (a basis penalty), and delta is the smoothing
## parameter selected by cross-validation inside the algorithm.  The fitted
## values are therefore a *linear smoother* of the working responses with hat
## matrix H = F (A + delta P)^{-1} t(F) Sigma^{-1}, and the effective number
## of parameters of that smoother is the familiar trace
##        edf = tr(H) = tr( (A + delta P)^{-1} A ) ,
## which lies between 0 (infinite penalty) and the nominal number of active
## coefficients (no penalty).  This is the standard penalized-regression /
## smoothing-spline definition of EDF (e.g. Hastie & Tibshirani; Wood, mgcv).
##
## Log-likelihood for the modified AIC.  At convergence the EM sets every
## residual covariance to the residual outer-product average of the *same*
## fitted coefficients, i.e.
##        Sigma_hat = (1/N) sum_i E[ r_i r_i' ] ,
## so the Gaussian quadratic form collapses to its profile value
##        sum_i E[ r_i' Sigma_hat^{-1} r_i ] = N * dim(r) .
## The (Monte-Carlo) complete-data Gaussian log-likelihood at the returned
## parameters is therefore available in closed form from the final variance
## parameters alone:
##        logLik = -1/2 sum_blocks N_b ( d_b log(2 pi) + log|Sigma_b| + d_b ) ,
## summed over the per-indicator latent-residual blocks (Sigma.fac, d=n.b),
## the scalar measurement-error blocks (Sigma.error, d=1, N_b = #obs of that
## indicator), and the per-latent structural blocks (Sigma.sem, d=n.b).
##
## The modified information criteria then pair this log-likelihood with the
## effective parameter count:
##        AIC.edf = -2 logLik + 2 * edf.total ,
##        BIC.edf = -2 logLik + log(N) * edf.total ,
## where edf.total = edf.coef (smoother traces) + edf.var (variance params).
## A classical AIC using the nominal coefficient count is also returned for
## reference.  The block is fully guarded: any failure leaves the estimate
## untouched and returns NULL.
## =====================================================================
.fsem_effective_df <- function(model, fit, n.b, data) {

  ## --- safe getter for the final variance parameters -----------------
  pe <- fit$result$params.estimated
  getpar <- function(name) {
    if (!is.null(pe$totparam) && !is.null(pe$totparam[[name]])) return(pe$totparam[[name]])
    pe[[name]]
  }
  sigma.fac   <- getpar("sigma.fac")
  sigma.error <- getpar("sigma.error")
  sigma.sem   <- getpar("sigma.sem")

  hist.fac <- fit$hist$hist.fac
  hist.sem <- fit$hist$hist.sem

  no.fac <- length(model$mod$factorModel)
  latents <- model$var$latents
  no.lat <- length(latents)

  symm    <- function(M) (M + t(M)) / 2
  log.det <- function(M) as.numeric(determinant(symm(as.matrix(M)), logarithm = TRUE)$modulus)

  ## --- 1. EDF of the factor-loading smoothers ------------------------
  edf.fac <- rep(NA_real_, no.fac)
  npar.fac <- rep(NA_real_, no.fac)
  for (j in seq_len(no.fac)) {
    hf <- hist.fac[[j]]
    if (is.null(hf) || is.null(hf$f) || is.null(hf$pen.fac)) next
    mont <- length(hf$f)
    Si <- solve(symm(as.matrix(sigma.fac[[j]])))
    p  <- ncol(hf$f[[1]][[1]])
    A  <- matrix(0, p, p)
    for (m in seq_len(mont)) {
      fm <- hf$f[[m]]
      for (k in seq_along(fm)) {
        Fk <- fm[[k]]
        A  <- A + crossprod(Fk, Si %*% Fk)
      }
    }
    A <- A / mont
    lam.sum <- A + hf$deltaa * as.matrix(hf$pen.fac)
    edf.fac[j]  <- sum(diag(solve(lam.sum, A)))
    npar.fac[j] <- p
  }

  ## --- 2. EDF of the structural-coefficient smoothers ----------------
  edf.sem  <- rep(0, no.lat)
  npar.sem <- rep(0, no.lat)
  names(edf.sem) <- names(npar.sem) <- latents
  for (lat in seq_len(no.lat)) {
    hs <- if (lat <= length(hist.sem)) hist.sem[[lat]] else NULL
    if (is.null(hs) || is.null(hs$g) || is.null(hs$pen.sem)) next
    mont <- length(hs$g)
    Si <- solve(symm(as.matrix(sigma.sem[[lat]])))
    g1 <- NULL
    for (m in seq_len(mont)) { for (i1 in seq_along(hs$g[[m]])) if (!is.null(hs$g[[m]][[i1]])) { g1 <- hs$g[[m]][[i1]]; break } ; if (!is.null(g1)) break }
    if (is.null(g1)) next
    p  <- ncol(g1)
    A  <- matrix(0, p, p)
    for (m in seq_len(mont)) {
      gm <- hs$g[[m]]
      for (i1 in seq_along(gm)) {
        Gk <- gm[[i1]]
        if (is.null(Gk)) next
        A <- A + crossprod(Gk, Si %*% Gk)
      }
    }
    A <- A / mont
    gam.sum <- A + hs$deltaa * as.matrix(hs$pen.sem)
    edf.sem[lat]  <- sum(diag(solve(gam.sum, A)))
    npar.sem[lat] <- p
  }

  ## --- 3. Variance/covariance parameter count (unpenalized) ----------
  cov.dim <- n.b * (n.b + 1) / 2
  edf.var <- no.fac * (cov.dim + 1) + no.lat * cov.dim

  ## --- 4. Sample sizes -----------------------------------------------
  N <- length(unique(data$.id))
  if (is.na(N) || N == 0) N <- length(hist.fac[[1]]$f[[1]])
  ## #observations contributing to each indicator's measurement error.
  ## `.ind` may store the indicator name (e.g. "z1") or its integer index.
  obs.count.ind <- vapply(seq_len(no.fac), function(j) {
    ind.name <- model$mod$factorModel[[j]]$indicator
    n.by.name <- sum(data$.ind == ind.name)
    if (n.by.name > 0) n.by.name else sum(data$.ind == j)
  }, numeric(1))

  ## --- 5. Closed-form (profile) complete-data Gaussian log-likelihood -
  l2pi <- log(2 * pi)
  loglik <- 0
  for (j in seq_len(no.fac)) {
    ld <- log.det(sigma.fac[[j]])
    loglik <- loglik - 0.5 * N * (n.b * l2pi + ld + n.b)
    s2 <- as.numeric(sigma.error[[j]])
    Ntot <- obs.count.ind[j]
    loglik <- loglik - 0.5 * Ntot * (l2pi + log(s2) + 1)
  }
  for (lat in seq_len(no.lat)) {
    if (is.null(sigma.sem[[lat]])) next
    ld <- log.det(sigma.sem[[lat]])
    loglik <- loglik - 0.5 * N * (n.b * l2pi + ld + n.b)
  }

  edf.coef  <- sum(edf.fac, na.rm = TRUE) + sum(edf.sem, na.rm = TRUE)
  npar.coef <- sum(npar.fac, na.rm = TRUE) + sum(npar.sem, na.rm = TRUE)
  edf.total  <- edf.coef + edf.var
  npar.total <- npar.coef + edf.var

  list(
    edf.fac    = edf.fac,
    edf.sem    = edf.sem,
    edf.coef   = edf.coef,
    edf.var    = edf.var,
    edf.total  = edf.total,
    npar.coef  = npar.coef,
    npar.total = npar.total,
    N          = N,
    loglik     = loglik,
    AIC        = -2 * loglik + 2 * npar.total,
    AIC.edf    = -2 * loglik + 2 * edf.total,
    BIC.edf    = -2 * loglik + log(N) * edf.total
  )
}

em.estimation.optimised.revised <- function(model, data, x.data = NULL,
                                            initial.parameter, n.b,
                                            n.em = 100, n.monte = 100,
                                            s.p = c("min", "sd"),
                                            range.min = NULL, range.max = NULL,
                                            design = c("regular", "irregular",
                                                       "regular.truncated", "regular.missing"),
                                            plot.progress = FALSE,
                                            save.every = NULL,
                                            gamma.init.sd = 0.3,
                                            seed = NULL) {

  init.perturbed <- initial.parameter
  if (!is.null(init.perturbed$gamma.param) && gamma.init.sd > 0) {
    if (!is.null(seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        old.seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit(assign(".Random.seed", old.seed, envir = .GlobalEnv), add = TRUE)
      }
      set.seed(seed)
    }
    for (k in seq_along(init.perturbed$gamma.param)) {
      g <- init.perturbed$gamma.param[[k]]
      if (length(g) > 0) {
        init.perturbed$gamma.param[[k]] <- g + rnorm(length(g), 0, gamma.init.sd)
      }
    }
  }

  fit <- em.estimation.optimised(model = model,
                          data = data,
                          x.data = x.data,
                          initial.parameter = init.perturbed,
                          n.b = n.b,
                          n.em = n.em,
                          n.monte = n.monte,
                          s.p = s.p,
                          range.min = range.min,
                          range.max = range.max,
                          design = design,
                          plot.progress = plot.progress,
                          save.every = save.every)

  ## ---- EDF + EDF-based modified AIC (separate block, never fatal) ----
  fit$edf <- tryCatch(
    .fsem_effective_df(model = model, fit = fit, n.b = n.b, data = data),
    error = function(e) {
      message("EDF computation skipped: ", conditionMessage(e))
      NULL
    }
  )

  fit
}
