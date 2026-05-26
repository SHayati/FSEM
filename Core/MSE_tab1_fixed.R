################################################################################
## MSE_tab1_fixed.R
##
## Compute Monte-Carlo MSE tables (factor part and SEM part) for the simulation
## produced by main_simulation_2_parallel.R.
##
## This is a corrected rewrite of Core/MSE_tab1.R. The original script had
## several issues that caused NaN / wrong MSE values:
##   * Stored norm entry was `sigma.error`, but the table-building loop read
##     `sigma.err` -> all sigma MSEs were NaN.
##   * In the inner factor / regression lookup, the loop variable `k` was used
##     instead of the model index `kk` when accessing $factor / $covariate.
##   * `d` (number of eigen components) leaked from one loop into the next,
##     so the sigma.sem.eigvec / eigval MSEs used the wrong dimension.
##   * The tabl.sem extraction loop iterated `1:length(mse[[l]][[nam.sem[j]]])`
##     (which includes intercept/eig.* columns) and tried to index `namm[m]`
##     out of bounds -> NA entries.
##   * `results` had to be in memory; it is now auto-loaded from
##     outputs/result_*_N_*_M_*.rds when missing.
##
## Output (in the calling environment):
##   tabl.fac : data.frame of factor-model MSEs (one row per N,M combination)
##   tabl.sem : data.frame of SEM-model MSEs   (one row per N,M combination)
##   mse, norm: intermediate per-run / per-cell objects
################################################################################

suppressPackageStartupMessages({
  library(Matrix)
})

## ---------------------------------------------------------------------------
## 1. Discover the result files (streamed; not held in memory together)
## ---------------------------------------------------------------------------
## The script supports two modes:
##   (a) `results` already in memory -> iterate over the list.
##   (b) Otherwise -> stream result_*.rds files from outputs/ one at a time.
## Mode (b) is required because each result file is hundreds of MB.
if (exists("results") && length(results) > 0) {
  use_files <- FALSE
  n <- length(results)
  message("Using in-memory `results` (", n, " runs).")
  get_run <- function(i) results[[i]]
  bootstrap <- results[[1]]
} else {
  use_files <- TRUE
  files <- list.files("outputs", pattern = "^result_\\d+_N_\\d+_M_\\d+\\.rds$",
                      full.names = TRUE)
  if (length(files) == 0)
    stop("No result_*.rds files found in 'outputs/'.")
  ## Optional sub-sampling for fast verification:
  ##   Sys.setenv(MSE_MAX_PER_CELL = "4")
  ##   source("Core/MSE_tab1_fixed.R")
  max_per_cell <- suppressWarnings(as.integer(Sys.getenv("MSE_MAX_PER_CELL", "")))
  if (!is.na(max_per_cell) && max_per_cell > 0) {
    key <- sub(".*_N_(\\d+)_M_(\\d+)\\.rds$", "\\1_\\2", basename(files))
    files <- unlist(lapply(split(files, key), function(g) head(g, max_per_cell)))
    message("Subsampling to ", max_per_cell, " runs per (N,M) cell -> ",
            length(files), " files.")
  }
  n <- length(files)
  message("Streaming ", n, " result files from outputs/ ...")
  get_run <- function(i) readRDS(files[i])
  bootstrap <- get_run(1)
}

## ---------------------------------------------------------------------------
## 2. Model bookkeeping (assumed identical across runs)
## ---------------------------------------------------------------------------
model.fit0  <- bootstrap$model.fit
n.coef.fac  <- length(model.fit0$var$indicators)
n.coef.sem  <- length(model.fit0$var$factors)
nam.fac     <- model.fit0$var$indicators
nam.sem     <- model.fit0$var$latents

## ---------------------------------------------------------------------------
## 3. Helper: sign-align estimated eigenvectors with the original ones
## ---------------------------------------------------------------------------
align_eigvecs <- function(est_list, org_list) {
  for (k in seq_along(est_list)) {
    Aest <- est_list[[k]]
    Aorg <- org_list[[k]]
    n.col <- ncol(Aest)
    n.row <- nrow(Aest)
    for (c in seq_len(n.col)) {
      s <- sum(Aest[, c] * Aorg[, c]) / n.row
      if (is.finite(s) && s < 0) Aest[, c] <- -Aest[, c]
    }
    est_list[[k]] <- Aest
  }
  est_list
}

## ---------------------------------------------------------------------------
## 4. Helper: find the effect type for an (indicator, factor) or
##    (response, covariate) pair in the fitted model
## ---------------------------------------------------------------------------
get_effect_fac <- function(model.fit, indicator, factor) {
  for (kk in seq_along(model.fit$mod$factorModel)) {
    fm <- model.fit$mod$factorModel[[kk]]
    if (identical(fm$indicator, indicator)) {
      idx <- which(fm$factor == factor)
      if (length(idx) > 0) return(fm$effect[idx[1]])
    }
  }
  NA_character_
}

get_effect_sem <- function(model.fit, response, covariate) {
  for (kk in seq_along(model.fit$mod$regression)) {
    rg <- model.fit$mod$regression[[kk]]
    if (identical(rg$response, response)) {
      idx <- which(rg$covariate == covariate)
      if (length(idx) > 0) return(rg$effect[idx[1]])
    }
  }
  NA_character_
}

## ---------------------------------------------------------------------------
## 5. Per-run squared distance (norm) between estimated and original params
## ---------------------------------------------------------------------------
norm <- vector("list", n)
names(norm) <- paste0("res", seq_len(n))

for (i in seq_len(n)) {
  if (use_files) message("  [", i, "/", n, "] ", basename(files[i]))
  res <- get_run(i)

  ## sample / time-grid sizes
  N <- length(unique(res$simulation$data$.id))
  M <- max(sapply(seq_len(N), function(j) {
    max(sapply(seq_along(res$model.fit$mod$factorModel), function(jj) {
      length(res$simulation$data$.t[res$simulation$data$.id == j &
                                    res$simulation$data$.ind == jj])
    }))
  }))

  ## sign-align eigenvectors before computing L2 distances
  res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std <-
    align_eigvecs(res$estimation$result$params.estimated.eval$sigma.fac.eigvec.std,
                  res$simulation$params.org.eval$sigma.fac.eigvec.std)
  res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std <-
    align_eigvecs(res$estimation$result$params.estimated.eval$sigma.sem.eigvec.std,
                  res$simulation$params.org.eval$sigma.sem.eigvec.std)

  est <- res$estimation$result$params.estimated.eval
  org <- res$simulation$params.org.eval

  cur <- list(N = N, M = M,
              coef.fac = list(), coef.sem = list(),
              intercept = list(),
              sigma.fac.eigvec = list(), sigma.sem.eigvec = list(),
              sigma.fac.eigval = list(), sigma.sem.eigval = list(),
              sigma.error = list())

  ## ---- factor loading coefficients -----------------------------------------
  for (j in seq_len(n.coef.fac)) {
    nm_j <- nam.fac[j]
    namm <- names(est$coef.fac.std[[j]])
    if (length(namm) > 0) {
      for (k in seq_along(namm)) {
        eff <- get_effect_fac(res$model.fit, nm_j, namm[k])
        diff_sq <- (est$coef.fac.std[[j]][[k]] - org$coef.fac.std[[j]][[k]])^2
        denom <- if (identical(eff, "historical")) 200^2 else 200
        cur$coef.fac[[nm_j]][[namm[k]]] <- sum(diff_sq) / denom
      }
    }
  }

  ## ---- structural regression coefficients ----------------------------------
  for (j in seq_len(n.coef.sem)) {
    nm_j <- nam.sem[j]
    namm <- names(est$coef.sem.std[[nm_j]])
    if (length(namm) > 0) {
      for (k in seq_along(namm)) {
        eff <- get_effect_sem(res$model.fit, nm_j, namm[k])
        diff_sq <- (est$coef.sem.std[[nm_j]][[k]] -
                    org$coef.sem.std[[nm_j]][[k]])^2
        denom <- if (identical(eff, "historical")) 200^2 else 200
        cur$coef.sem[[nm_j]][[namm[k]]] <- sum(diff_sq) / denom
      }
    }
  }

  ## ---- intercepts ----------------------------------------------------------
  for (j in seq_len(n.coef.fac)) {
    nm_j <- nam.fac[j]
    cur$intercept[[nm_j]] <-
      sum((est$intercept.std[[j]] - org$intercept.std[[j]])^2) / 200
  }

  ## ---- factor-side covariance eigen-components -----------------------------
  for (j in seq_len(n.coef.fac)) {
    nm_j <- nam.fac[j]
    d.vec <- ncol(est$sigma.fac.eigvec.std[[nm_j]])
    cur$sigma.fac.eigvec[[nm_j]] <- numeric(d.vec)
    cur$sigma.fac.eigval[[nm_j]] <- numeric(d.vec)
    for (k in seq_len(d.vec)) {
      cur$sigma.fac.eigvec[[nm_j]][k] <-
        sum((est$sigma.fac.eigvec.std[[j]][, k] -
             org$sigma.fac.eigvec.std[[j]][, k])^2) / 200
      cur$sigma.fac.eigval[[nm_j]][k] <-
        (est$sigma.fac.eigval.std[[j]][k] -
         org$sigma.fac.eigval.std[[j]][k])^2
    }
  }

  ## ---- structural covariance eigen-components ------------------------------
  for (j in seq_len(n.coef.sem)) {
    nm_j <- nam.sem[j]
    d.vec <- ncol(est$sigma.sem.eigvec.std[[nm_j]])
    cur$sigma.sem.eigvec[[nm_j]] <- numeric(d.vec)
    cur$sigma.sem.eigval[[nm_j]] <- numeric(d.vec)
    for (k in seq_len(d.vec)) {
      cur$sigma.sem.eigvec[[nm_j]][k] <-
        sum((est$sigma.sem.eigvec.std[[j]][, k] -
             org$sigma.sem.eigvec.std[[j]][, k])^2) / 200
      cur$sigma.sem.eigval[[nm_j]][k] <-
        (est$sigma.sem.eigval.std[[j]][k] -
         org$sigma.sem.eigval.std[[j]][k])^2
    }
  }

  ## ---- measurement error variances -----------------------------------------
  for (j in seq_len(n.coef.fac)) {
    nm_j <- nam.fac[j]
    cur$sigma.error[[nm_j]] <-
      (res$estimation$result$params.estimated$totparam$sigma.error[[j]] -
       org$sigma.error[[j]])^2
  }

  norm[[i]] <- cur
  rm(res); gc(verbose = FALSE)
}

## ---------------------------------------------------------------------------
## 6. Average over Monte-Carlo runs within each (N, M) cell
## ---------------------------------------------------------------------------
comb <- expand.grid(N = c(50, 100), M = c(10, 20))
comb <- comb[order(comb$N, comb$M), ]
rownames(comb) <- NULL

est0 <- bootstrap$estimation$result$params.estimated.eval
rm(bootstrap); gc(verbose = FALSE)

mean_or_na <- function(x) if (length(x) == 0) NA_real_ else mean(x, na.rm = TRUE)

mse <- vector("list", nrow(comb))
for (l in seq_len(nrow(comb))) {
  Nl <- comb$N[l]; Ml <- comb$M[l]
  sel <- which(vapply(norm, function(nn) nn$N == Nl && nn$M == Ml, logical(1)))
  cell <- list(N = Nl, M = Ml, numb.runs = length(sel))

  ## factor side
  for (j in seq_len(n.coef.fac)) {
    nm_j <- nam.fac[j]
    namm <- names(est0$coef.fac.std[[j]])
    sub  <- list()
    if (length(namm) > 0) {
      for (k in seq_along(namm)) {
        vals <- vapply(sel, function(i) {
          v <- norm[[i]]$coef.fac[[nm_j]][[namm[k]]]
          if (is.null(v)) NA_real_ else v
        }, numeric(1))
        sub[[namm[k]]] <- mean_or_na(vals)
      }
    }
    sub$intercept <- mean_or_na(vapply(sel, function(i)
      norm[[i]]$intercept[[nm_j]], numeric(1)))
    for (k in 1:2) {
      sub[[paste0("eig.fuc.", k)]] <- mean_or_na(vapply(sel, function(i)
        norm[[i]]$sigma.fac.eigvec[[nm_j]][k], numeric(1)))
      sub[[paste0("eig.val.", k)]] <- mean_or_na(vapply(sel, function(i)
        norm[[i]]$sigma.fac.eigval[[nm_j]][k], numeric(1)))
    }
    sub$sig.err <- mean_or_na(vapply(sel, function(i)
      norm[[i]]$sigma.error[[nm_j]], numeric(1)))
    cell[[nm_j]] <- sub
  }

  ## SEM side
  for (j in seq_len(n.coef.sem)) {
    nm_j <- nam.sem[j]
    namm <- names(est0$coef.sem.std[[nm_j]])
    sub  <- list()
    if (length(namm) > 0) {
      for (k in seq_along(namm)) {
        vals <- vapply(sel, function(i) {
          v <- norm[[i]]$coef.sem[[nm_j]][[namm[k]]]
          if (is.null(v)) NA_real_ else v
        }, numeric(1))
        sub[[namm[k]]] <- mean_or_na(vals)
      }
    }
    for (k in 1:2) {
      sub[[paste0("eig.fuc.", k)]] <- mean_or_na(vapply(sel, function(i)
        norm[[i]]$sigma.sem.eigvec[[nm_j]][k], numeric(1)))
      sub[[paste0("eig.val.", k)]] <- mean_or_na(vapply(sel, function(i)
        norm[[i]]$sigma.sem.eigval[[nm_j]][k], numeric(1)))
    }
    cell[[nm_j]] <- sub
  }

  mse[[l]] <- cell
}

## ---------------------------------------------------------------------------
## 7. Build the two MSE tables
## ---------------------------------------------------------------------------
tabl.fac <- data.frame(N = comb$N, M = comb$M)
for (l in seq_len(nrow(comb))) {
  for (j in seq_len(n.coef.fac)) {
    nm_j <- nam.fac[j]
    tabl.fac[[paste0("beta",   j)]][l] <- mse[[l]][[nm_j]]$intercept
    namm <- names(est0$coef.fac.std[[j]])
    if (length(namm) > 0) {
      for (m in seq_along(namm))
        tabl.fac[[paste0("lambda", j, m)]][l] <- mse[[l]][[nm_j]][[namm[m]]]
    }
    tabl.fac[[paste0("phi",   j, 1)]][l] <- mse[[l]][[nm_j]]$eig.fuc.1
    tabl.fac[[paste0("phi",   j, 2)]][l] <- mse[[l]][[nm_j]]$eig.fuc.2
    tabl.fac[[paste0("nu",    j, 1)]][l] <- mse[[l]][[nm_j]]$eig.val.1
    tabl.fac[[paste0("nu",    j, 2)]][l] <- mse[[l]][[nm_j]]$eig.val.2
    tabl.fac[[paste0("sigma", j)]][l]    <- mse[[l]][[nm_j]]$sig.err
  }
}

tabl.sem <- data.frame(N = comb$N, M = comb$M)
for (l in seq_len(nrow(comb))) {
  for (j in seq_len(n.coef.sem)) {
    nm_j <- nam.sem[j]
    namm <- names(est0$coef.sem.std[[nm_j]])
    if (length(namm) > 0) {
      for (m in seq_along(namm))
        tabl.sem[[paste0("gamma", j, m)]][l] <- mse[[l]][[nm_j]][[namm[m]]]
    }
    tabl.sem[[paste0("psi", j, 1)]][l] <- mse[[l]][[nm_j]]$eig.fuc.1
    tabl.sem[[paste0("psi", j, 2)]][l] <- mse[[l]][[nm_j]]$eig.fuc.2
    tabl.sem[[paste0("mu",  j, 1)]][l] <- mse[[l]][[nm_j]]$eig.val.1
    tabl.sem[[paste0("mu",  j, 2)]][l] <- mse[[l]][[nm_j]]$eig.val.2
  }
}

## ---------------------------------------------------------------------------
## 8. Persist tables to outputs/
## ---------------------------------------------------------------------------
if (!dir.exists("outputs")) dir.create("outputs")
write.csv(tabl.fac, file.path("outputs", "MSE_tab1_fac.csv"), row.names = FALSE)
write.csv(tabl.sem, file.path("outputs", "MSE_tab1_sem.csv"), row.names = FALSE)

cat("\n=== MSE table (factor part) ===\n");      print(tabl.fac)
cat("\n=== MSE table (structural part) ===\n");  print(tabl.sem)
