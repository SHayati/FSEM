# ============================================================================
# Compare estimation quality for simulation_3 (dependent latents) vs
# simulation_2_parallel (uncorrelated latents).
#
# Notes:
# - parameters are not identifiable in raw form; we only compare the
#   STANDARDIZED versions (`*.std`) which ARE identifiable.
# - Sign of estimated eigenfunctions is flipped if it points opposite the
#   truth (same convention as Core/MSE_tab1.R).
# - L2 norm uses the same Riemann grid normalizer (200) as MSE_tab1.R;
#   historical 2D kernels use 200^2.
# ============================================================================

suppressPackageStartupMessages({
  library(Matrix)
})

n.grid <- 200  # Riemann grid used inside simulation()

# -------- helpers ------------------------------------------------------------
align_sign <- function(est, tru) {
  d <- dim(est)
  for (col in 1:d[2]) {
    s <- mean(est[, col] * tru[, col])
    if (s < 0) est[, col] <- -est[, col]
  }
  est
}

l2_vec <- function(a, b) mean((a - b)^2)                  # 1/n.grid * sum
l2_mat <- function(A, B) mean((A - B)^2)                  # 1/n.grid^2 * sum

# ---- per-result MSEs for a single (sim, est) pair --------------------------
mse_one <- function(sim, est_eval, model_fit, est_full = NULL) {
  ind <- model_fit$var$indicators
  lat <- model_fit$var$latents

  out <- list()

  # -- factor-model coefficients (lambda) & intercepts (beta) ---------------
  for (j in seq_along(ind)) {
    nm_j <- ind[j]
    out$lambda[[nm_j]] <- list()

    fm_idx <- which(vapply(model_fit$mod$factorModel,
                           function(f) f$indicator == nm_j, logical(1)))
    fm <- model_fit$mod$factorModel[[fm_idx]]

    nms_k <- names(est_eval$coef.fac.std[[j]])
    for (k in seq_along(nms_k)) {
      mm <- which(fm$factor == nms_k[k])
      eff <- fm$effect[mm]
      e <- est_eval$coef.fac.std[[j]][[k]]
      t <- sim$params.org.eval$coef.fac.std[[j]][[k]]
      out$lambda[[nm_j]][[nms_k[k]]] <-
        if (eff == "historical") l2_mat(e, t) else l2_vec(e, t)
    }

    out$beta[[nm_j]] <- l2_vec(est_eval$intercept.std[[j]],
                               sim$params.org.eval$intercept.std[[j]])

    # eigen-functions / values of factor error covariance
    Ev <- align_sign(est_eval$sigma.fac.eigvec.std[[j]],
                     sim$params.org.eval$sigma.fac.eigvec.std[[j]])
    out$phi[[nm_j]] <- c(l2_vec(Ev[, 1], sim$params.org.eval$sigma.fac.eigvec.std[[j]][, 1]),
                         l2_vec(Ev[, 2], sim$params.org.eval$sigma.fac.eigvec.std[[j]][, 2]))
    out$nu[[nm_j]]  <- (est_eval$sigma.fac.eigval.std[[j]][1:2] -
                          sim$params.org.eval$sigma.fac.eigval.std[[j]][1:2])^2

    if (!is.null(est_full)) {
      out$sig_err[[nm_j]] <-
        (est_full$result$params.estimated$totparam$sigma.error[[j]] -
           sim$params.org.eval$sigma.error[[j]])^2
    }
  }

  # -- structural (SEM) coefficients (gamma) --------------------------------
  for (j in seq_along(lat)) {
    nm_j <- lat[j]
    nms_k <- names(est_eval$coef.sem.std[[nm_j]])
    if (length(nms_k) > 0) {
      reg_idx <- which(vapply(model_fit$mod$regression,
                              function(r) r$response == nm_j, logical(1)))
      reg <- model_fit$mod$regression[[reg_idx]]
      out$gamma[[nm_j]] <- list()
      for (k in seq_along(nms_k)) {
        mm <- which(reg$covariate == nms_k[k])
        eff <- reg$effect[mm]
        e <- est_eval$coef.sem.std[[nm_j]][[k]]
        t <- sim$params.org.eval$coef.sem.std[[nm_j]][[k]]
        out$gamma[[nm_j]][[nms_k[k]]] <-
          if (eff == "historical") l2_mat(e, t) else l2_vec(e, t)
      }
    }

    # eigen-functions / values of structural error covariance
    Ev <- align_sign(est_eval$sigma.sem.eigvec.std[[j]],
                     sim$params.org.eval$sigma.sem.eigvec.std[[j]])
    out$psi[[nm_j]] <- c(l2_vec(Ev[, 1], sim$params.org.eval$sigma.sem.eigvec.std[[j]][, 1]),
                         l2_vec(Ev[, 2], sim$params.org.eval$sigma.sem.eigvec.std[[j]][, 2]))
    out$mu[[nm_j]]  <- (est_eval$sigma.sem.eigval.std[[j]][1:2] -
                          sim$params.org.eval$sigma.sem.eigval.std[[j]][1:2])^2
  }
  out
}

mse_one_wrap <- function(res) {
  est <- res$estimation
  mse_one(res$simulation, est$result$params.estimated.eval, res$model.fit, est)
}

# ============================================================================
# 1) Single run from main_simulation_3.R
# ============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("MAIN_SIMULATION_3: single run (dependent factors, eta1 ~ eta2 concurrent)\n")
cat(strrep("=", 70), "\n", sep = "")

sr <- readRDS("single_run_results.rds")

# Both fitted models use model.sim's variable layout
res_w <- list(simulation = sr$sim.data,
              estimation = sr$params.estimated.well.specified,
              model.fit  = sr$sim.data$params$model)

# need a model.fit list with $var and $mod. params$model is just $factorModel,$regression.
# Reconstruct by re-evaluating fsem() expressions.
source("Core/FSEM.R")
source("Core/simualtion_utilities.R")

mk_models <- function() {
  model.sim <- fsem(eta1 ~~ z1 + z2,                effectType = "concurrent") %+%
               fsem(eta1 ~~ z3,                     effectType = "historical") %+%
               fsem(eta2 ~~ z4 + z5 + z6,           effectType = "concurrent") %+%
               fsem(eta1 ~ -1 + eta2,               effectType = "concurrent",
                    latent.covariate = "eta2", scalar.covariate = FALSE) %+%
               fsem(eta2 ~ -1)
  model.fit.m <- fsem(eta1 ~~ z1 + z2,              effectType = "concurrent") %+%
                 fsem(eta1 ~~ z3,                   effectType = "historical") %+%
                 fsem(eta2 ~~ z4 + z5 + z6,         effectType = "concurrent") %+%
                 fsem(eta1 ~ -1,                    effectType = "concurrent",
                      latent.covariate = "eta2", scalar.covariate = FALSE) %+%
                 fsem(eta2 ~ -1)
  list(w = model.sim, m = model.fit.m)
}
mods <- mk_models()

m_w <- mse_one(sr$sim.data,
               sr$params.estimated.well.specified$result$params.estimated.eval,
               mods$w, sr$params.estimated.well.specified)
m_m <- mse_one(sr$sim.data,
               sr$params.estimated.misspecified$result$params.estimated.eval,
               mods$m, sr$params.estimated.misspecified)

print_one <- function(tag, m) {
  cat("\n----", tag, "----\n")
  cat("lambda (factor loadings) MSE per indicator/factor:\n")
  for (ind in names(m$lambda))
    for (f in names(m$lambda[[ind]]))
      cat(sprintf("  %s ~ %s : %.5f\n", ind, f, m$lambda[[ind]][[f]]))
  cat("beta (intercepts) MSE per indicator:\n")
  for (ind in names(m$beta)) cat(sprintf("  %s : %.5f\n", ind, m$beta[[ind]]))
  cat("phi (factor-error eigfuns) MSE per indicator [ef1, ef2]:\n")
  for (ind in names(m$phi))
    cat(sprintf("  %s : %.5f, %.5f\n", ind, m$phi[[ind]][1], m$phi[[ind]][2]))
  cat("nu (factor-error eigvals) MSE per indicator [ev1, ev2]:\n")
  for (ind in names(m$nu))
    cat(sprintf("  %s : %.5f, %.5f\n", ind, m$nu[[ind]][1], m$nu[[ind]][2]))
  cat("sigma_err MSE per indicator:\n")
  for (ind in names(m$sig_err))
    cat(sprintf("  %s : %.5f\n", ind, m$sig_err[[ind]]))
  cat("gamma (structural coefs) MSE per latent/covariate:\n")
  if (!is.null(m$gamma))
    for (lt in names(m$gamma))
      for (cv in names(m$gamma[[lt]]))
        cat(sprintf("  %s ~ %s : %.5f\n", lt, cv, m$gamma[[lt]][[cv]]))
  cat("psi (struct-error eigfuns) MSE per latent [ef1, ef2]:\n")
  for (lt in names(m$psi))
    cat(sprintf("  %s : %.5f, %.5f\n", lt, m$psi[[lt]][1], m$psi[[lt]][2]))
  cat("mu (struct-error eigvals) MSE per latent [ev1, ev2]:\n")
  for (lt in names(m$mu))
    cat(sprintf("  %s : %.5f, %.5f\n", lt, m$mu[[lt]][1], m$mu[[lt]][2]))
}

print_one("WELL-SPECIFIED FIT",  m_w)
print_one("MISSPECIFIED  FIT (omits eta1<-eta2)", m_m)

# ============================================================================
# 2) Reference: monte-carlo results from main_simulation_2_parallel.R
# (uncorrelated factors; 4 (N,M) combos with multiple replicates each)
# ============================================================================
cat("\n\n", strrep("=", 70), "\n", sep = "")
cat("MAIN_SIMULATION_2_PARALLEL: monte-carlo (uncorrelated factors)\n")
cat(strrep("=", 70), "\n", sep = "")

files <- list.files("outputs", pattern = "^result_.*\\.rds$", full.names = TRUE)
cat("Found", length(files), "result files in outputs/.\n")

# group by (N,M)
meta <- do.call(rbind, lapply(files, function(f) {
  m <- regmatches(basename(f),
                  regexec("result_\\d+_N_(\\d+)_M_(\\d+)\\.rds", basename(f)))[[1]]
  data.frame(file = f, N = as.integer(m[2]), M = as.integer(m[3]),
             stringsAsFactors = FALSE)
}))

agg <- list()
for (k in seq_len(nrow(meta))) {
  r <- tryCatch(readRDS(meta$file[k]), error = function(e) NULL)
  if (is.null(r)) next
  if (is.null(r$model.fit) || is.null(r$estimation$result$params.estimated.eval))
    next
  mm <- tryCatch(mse_one_wrap(r), error = function(e) NULL)
  if (is.null(mm)) next
  key <- sprintf("N=%d,M=%d", meta$N[k], meta$M[k])
  agg[[key]] <- c(agg[[key]], list(mm))
}

summ_combo <- function(lst) {
  # average all numeric leaves across replicates
  collect <- function(x) {
    if (is.numeric(x)) return(list(x))
    if (is.list(x))    return(unlist(lapply(x, collect), recursive = FALSE))
    list()
  }
  flat <- lapply(lst, function(m) {
    keys <- c()
    vals <- c()
    walk <- function(node, path) {
      if (is.numeric(node)) {
        for (i in seq_along(node)) {
          keys <<- c(keys, if (length(node) > 1)
                      paste0(path, "[", i, "]") else path)
          vals <<- c(vals, node[i])
        }
        return()
      }
      if (is.list(node)) {
        nm <- names(node); if (is.null(nm)) nm <- as.character(seq_along(node))
        for (i in seq_along(node)) walk(node[[i]], paste0(path, "/", nm[i]))
      }
    }
    walk(m, "")
    setNames(vals, keys)
  })
  all_keys <- Reduce(union, lapply(flat, names))
  mat <- sapply(flat, function(v) v[all_keys])
  rownames(mat) <- all_keys
  rowMeans(mat, na.rm = TRUE)
}

cat("\nAveraged per-element MSE for each (N,M) combination:\n")
for (key in names(agg)) {
  cat("\n----", key, "  (replicates=", length(agg[[key]]), ") ----\n", sep = "")
  s <- summ_combo(agg[[key]])
  # show top-line summary by group
  groups <- list(
    lambda  = s[grep("^/lambda/", names(s))],
    beta    = s[grep("^/beta/",   names(s))],
    phi     = s[grep("^/phi/",    names(s))],
    nu      = s[grep("^/nu/",     names(s))],
    sig_err = s[grep("^/sig_err/",names(s))],
    gamma   = s[grep("^/gamma/",  names(s))],
    psi     = s[grep("^/psi/",    names(s))],
    mu      = s[grep("^/mu/",     names(s))]
  )
  for (g in names(groups))
    if (length(groups[[g]]))
      cat(sprintf("  %-7s : mean MSE = %.5f  (n elems=%d)\n",
                  g, mean(groups[[g]]), length(groups[[g]])))
}

# ============================================================================
# 3) Side-by-side comparison summary
# ============================================================================
cat("\n\n", strrep("=", 70), "\n", sep = "")
cat("SUMMARY: simulation_3 (single run) vs simulation_2 (avg MC, N=50,M=20)\n")
cat(strrep("=", 70), "\n", sep = "")

flatten_one <- function(m) {
  keys <- c(); vals <- c()
  walk <- function(node, path) {
    if (is.numeric(node)) {
      for (i in seq_along(node)) {
        keys <<- c(keys, if (length(node) > 1)
                    paste0(path, "[", i, "]") else path)
        vals <<- c(vals, node[i])
      }
      return()
    }
    if (is.list(node)) {
      nm <- names(node); if (is.null(nm)) nm <- as.character(seq_along(node))
      for (i in seq_along(node)) walk(node[[i]], paste0(path, "/", nm[i]))
    }
  }
  walk(m, "")
  setNames(vals, keys)
}

f_w  <- flatten_one(m_w)
f_m  <- flatten_one(m_m)

# closest reference combo (single-run used N=50, M=20)
ref_key <- "N=50,M=20"
ref <- if (!is.null(agg[[ref_key]])) summ_combo(agg[[ref_key]]) else NULL

groups <- c("lambda","beta","phi","nu","sig_err","gamma","psi","mu")
cat(sprintf("%-8s | %12s | %12s | %12s\n",
            "param", "sim3-wellsp", "sim3-misspc", "sim2(N=50,M=20)"))
cat(strrep("-", 60), "\n", sep = "")
for (g in groups) {
  patt <- paste0("^/", g, "/")
  vw <- mean(f_w[grep(patt, names(f_w))])
  vm <- mean(f_m[grep(patt, names(f_m))])
  vr <- if (!is.null(ref)) mean(ref[grep(patt, names(ref))]) else NA
  cat(sprintf("%-8s | %12.5f | %12.5f | %12.5f\n", g, vw, vm, vr))
}

cat("\nDone.\n")
