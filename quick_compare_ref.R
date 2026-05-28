suppressPackageStartupMessages({library(Matrix)})

l2_vec <- function(a, b) mean((a - b)^2)
l2_mat <- function(A, B) mean((A - B)^2)

align_sign <- function(est, tru) {
  d <- dim(est)
  if (length(d) < 2) return(est)
  for (col in seq_len(d[2])) if (mean(est[, col] * tru[, col]) < 0) est[, col] <- -est[, col]
  est
}

safe_mean <- function(x) if (length(x)) mean(x, na.rm = TRUE) else NA_real_

score_one <- function(sim, est_eval, est_full = NULL) {
  out <- list(
    lambda = numeric(), beta = numeric(), phi = numeric(), nu = numeric(),
    sig_err = numeric(), gamma = numeric(), psi = numeric(), mu = numeric()
  )

  ind <- names(est_eval$coef.fac.std)
  for (nm in ind) {
    ef <- est_eval$coef.fac.std[[nm]]
    tf <- sim$params.org.eval$coef.fac.std[[nm]]
    if (!is.null(ef) && !is.null(tf)) {
      for (k in intersect(names(ef), names(tf))) {
        if (is.matrix(ef[[k]]) || is.matrix(tf[[k]])) out$lambda <- c(out$lambda, l2_mat(ef[[k]], tf[[k]]))
        else out$lambda <- c(out$lambda, l2_vec(ef[[k]], tf[[k]]))
      }
    }
    if (!is.null(est_eval$intercept.std[[nm]]) && !is.null(sim$params.org.eval$intercept.std[[nm]])) {
      out$beta <- c(out$beta, l2_vec(est_eval$intercept.std[[nm]], sim$params.org.eval$intercept.std[[nm]]))
    }
    if (!is.null(est_eval$sigma.fac.eigvec.std[[nm]]) && !is.null(sim$params.org.eval$sigma.fac.eigvec.std[[nm]])) {
      ev <- align_sign(est_eval$sigma.fac.eigvec.std[[nm]], sim$params.org.eval$sigma.fac.eigvec.std[[nm]])
      kmax <- min(2, ncol(ev), ncol(sim$params.org.eval$sigma.fac.eigvec.std[[nm]]))
      for (k in seq_len(kmax)) out$phi <- c(out$phi, l2_vec(ev[, k], sim$params.org.eval$sigma.fac.eigvec.std[[nm]][, k]))
    }
    if (!is.null(est_eval$sigma.fac.eigval.std[[nm]]) && !is.null(sim$params.org.eval$sigma.fac.eigval.std[[nm]])) {
      kmax <- min(2, length(est_eval$sigma.fac.eigval.std[[nm]]), length(sim$params.org.eval$sigma.fac.eigval.std[[nm]]))
      out$nu <- c(out$nu, (est_eval$sigma.fac.eigval.std[[nm]][1:kmax] - sim$params.org.eval$sigma.fac.eigval.std[[nm]][1:kmax])^2)
    }
  }

  lat <- names(est_eval$coef.sem.std)
  for (nm in lat) {
    es <- est_eval$coef.sem.std[[nm]]
    ts <- sim$params.org.eval$coef.sem.std[[nm]]
    if (!is.null(es) && !is.null(ts)) {
      for (k in intersect(names(es), names(ts))) {
        if (is.matrix(es[[k]]) || is.matrix(ts[[k]])) out$gamma <- c(out$gamma, l2_mat(es[[k]], ts[[k]]))
        else out$gamma <- c(out$gamma, l2_vec(es[[k]], ts[[k]]))
      }
    }
    if (!is.null(est_eval$sigma.sem.eigvec.std[[nm]]) && !is.null(sim$params.org.eval$sigma.sem.eigvec.std[[nm]])) {
      ev <- align_sign(est_eval$sigma.sem.eigvec.std[[nm]], sim$params.org.eval$sigma.sem.eigvec.std[[nm]])
      kmax <- min(2, ncol(ev), ncol(sim$params.org.eval$sigma.sem.eigvec.std[[nm]]))
      for (k in seq_len(kmax)) out$psi <- c(out$psi, l2_vec(ev[, k], sim$params.org.eval$sigma.sem.eigvec.std[[nm]][, k]))
    }
    if (!is.null(est_eval$sigma.sem.eigval.std[[nm]]) && !is.null(sim$params.org.eval$sigma.sem.eigval.std[[nm]])) {
      kmax <- min(2, length(est_eval$sigma.sem.eigval.std[[nm]]), length(sim$params.org.eval$sigma.sem.eigval.std[[nm]]))
      out$mu <- c(out$mu, (est_eval$sigma.sem.eigval.std[[nm]][1:kmax] - sim$params.org.eval$sigma.sem.eigval.std[[nm]][1:kmax])^2)
    }
  }

  if (!is.null(est_full) && !is.null(est_full$result$params.estimated$totparam$sigma.error) && !is.null(sim$params.org.eval$sigma.error)) {
    se <- est_full$result$params.estimated$totparam$sigma.error
    te <- sim$params.org.eval$sigma.error
    n <- min(length(se), length(te))
    if (n > 0) out$sig_err <- (unlist(se[1:n]) - unlist(te[1:n]))^2
  }

  sapply(out, safe_mean)
}

sr <- readRDS("single_run_results.rds")
well <- score_one(sr$sim.data, sr$params.estimated.well.specified$result$params.estimated.eval, sr$params.estimated.well.specified)
miss <- score_one(sr$sim.data, sr$params.estimated.misspecified$result$params.estimated.eval, sr$params.estimated.misspecified)

files <- list.files("outputs", pattern = "_N_50_M_20\\.rds$", full.names = TRUE)
cat("Found", length(files), "reference files for N=50, M=20\n")

scores <- list()
for (i in seq_along(files)) {
  r <- tryCatch(readRDS(files[i]), error = function(e) NULL)
  if (is.null(r) || is.null(r$simulation) || is.null(r$estimation$result$params.estimated.eval)) next
  sc <- tryCatch(score_one(r$simulation, r$estimation$result$params.estimated.eval, r$estimation), error = function(e) NULL)
  if (is.null(sc)) next
  scores[[length(scores) + 1]] <- sc
  if (i %% 20 == 0) cat("processed", i, "files\n")
}

ref <- colMeans(do.call(rbind, scores), na.rm = TRUE)

cat("\nSingle-run (dependent factors):\n")
print(rbind(sim3_well = well, sim3_miss = miss))
cat("\nReference MC average (uncorrelated factors), N=50 M=20:\n")
print(ref)
cat("\nRelative to reference (ratio, <1 better, >1 worse):\n")
print(rbind(well_over_ref = well / ref, miss_over_ref = miss / ref))
