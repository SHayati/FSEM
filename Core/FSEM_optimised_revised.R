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

  em.estimation.optimised(model = model,
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
}
