library(Matrix)

.fsem_align_eigenvectors <- function(estimated, truth) {
  if (length(estimated) == 0 || length(truth) == 0) {
    return(estimated)
  }

  for (index in seq_along(estimated)) {
    est_mat <- estimated[[index]]
    truth_mat <- truth[[index]]
    if (is.null(dim(est_mat)) || is.null(dim(truth_mat))) {
      next
    }

    n_cols <- min(ncol(est_mat), ncol(truth_mat))
    for (col_idx in seq_len(n_cols)) {
      inner_prod <- mean(est_mat[, col_idx] * truth_mat[, col_idx])
      if (!is.na(inner_prod) && inner_prod < 0) {
        estimated[[index]][, col_idx] <- -estimated[[index]][, col_idx]
      }
    }
  }

  estimated
}

.fsem_effect_type <- function(model_list, response_name, covariate_name, response_field, covariate_field) {
  for (entry in model_list) {
    if (!identical(entry[[response_field]], response_name)) {
      next
    }

    covariates <- entry[[covariate_field]]
    match_idx <- match(covariate_name, covariates)
    if (!is.na(match_idx)) {
      return(entry$effect[[match_idx]])
    }
  }

  NA_character_
}

.fsem_compute_single_mse <- function(dat_list, fit_key, estimation_key) {
  if (length(dat_list) == 0) {
    stop("results is empty")
  }

  first_res <- dat_list[[1]]
  model_fit <- first_res[[fit_key]]
  estimation <- first_res[[estimation_key]]

  nam.fac <- model_fit$var$indicators
  nam.sem <- model_fit$var$latents
  combo_df <- unique(do.call(rbind, lapply(dat_list, function(res) {
    data.frame(
      N = length(unique(res$simulation$data$.id)),
      M = max(sapply(unique(res$simulation$data$.id), function(id_value) {
        max(sapply(seq_along(model_fit$mod$factorModel), function(ind_idx) {
          length(res$simulation$data$.t[
            res$simulation$data$.id == id_value &
              res$simulation$data$.ind == ind_idx
          ])
        }))
      }))
    )
  })))
  combo_df <- combo_df[order(combo_df$N, combo_df$M), , drop = FALSE]
  rownames(combo_df) <- NULL

  norm_list <- vector("list", length(dat_list))

  for (res_idx in seq_along(dat_list)) {
    res <- dat_list[[res_idx]]
    est_eval <- res[[estimation_key]]$result$params.estimated.eval
    true_eval <- res$simulation$params.org.eval

    est_eval$sigma.fac.eigvec.std <- .fsem_align_eigenvectors(
      est_eval$sigma.fac.eigvec.std,
      true_eval$sigma.fac.eigvec.std
    )
    est_eval$sigma.sem.eigvec.std <- .fsem_align_eigenvectors(
      est_eval$sigma.sem.eigvec.std,
      true_eval$sigma.sem.eigvec.std
    )

    N <- length(unique(res$simulation$data$.id))
    M <- max(sapply(unique(res$simulation$data$.id), function(id_value) {
      max(sapply(seq_along(res[[fit_key]]$mod$factorModel), function(ind_idx) {
        length(res$simulation$data$.t[
          res$simulation$data$.id == id_value &
            res$simulation$data$.ind == ind_idx
        ])
      }))
    }))

    entry <- list(N = N, M = M)

    for (fac_idx in seq_along(nam.fac)) {
      indicator_name <- nam.fac[[fac_idx]]
      fac_names <- names(est_eval$coef.fac.std[[fac_idx]])
      entry$coef.fac[[indicator_name]] <- list()

      if (length(fac_names) > 0) {
        for (name_idx in seq_along(fac_names)) {
          factor_name <- fac_names[[name_idx]]
          effect_type <- .fsem_effect_type(
            res[[fit_key]]$mod$factorModel,
            indicator_name,
            factor_name,
            "indicator",
            "factor"
          )
          denom <- if (identical(effect_type, "historical")) 200^2 else 200
          entry$coef.fac[[indicator_name]][[factor_name]] <- sum(
            (est_eval$coef.fac.std[[fac_idx]][[name_idx]] -
               true_eval$coef.fac.std[[fac_idx]][[name_idx]])^2
          ) / denom
        }
      }

      entry$intercept[[indicator_name]] <- mean(
        (est_eval$intercept.std[[fac_idx]] - true_eval$intercept.std[[fac_idx]])^2
      )

      n_fac_cols <- ncol(est_eval$sigma.fac.eigvec.std[[fac_idx]])
      entry$sigma.fac.eigvec[[indicator_name]] <- vapply(seq_len(n_fac_cols), function(col_idx) {
        mean((est_eval$sigma.fac.eigvec.std[[fac_idx]][, col_idx] -
                true_eval$sigma.fac.eigvec.std[[fac_idx]][, col_idx])^2)
      }, numeric(1))
      entry$sigma.fac.eigval[[indicator_name]] <- vapply(seq_along(est_eval$sigma.fac.eigval.std[[fac_idx]]), function(col_idx) {
        (est_eval$sigma.fac.eigval.std[[fac_idx]][[col_idx]] -
           true_eval$sigma.fac.eigval.std[[fac_idx]][[col_idx]])^2
      }, numeric(1))
      entry$sigma.error[[indicator_name]] <- (
        res[[estimation_key]]$result$params.estimated$totparam$sigma.error[[fac_idx]] -
          res$simulation$params.org.eval$sigma.error[[fac_idx]]
      )^2
    }

    for (sem_idx in seq_along(nam.sem)) {
      latent_name <- nam.sem[[sem_idx]]
      sem_names <- names(est_eval$coef.sem.std[[latent_name]])
      entry$coef.sem[[latent_name]] <- list()

      if (length(sem_names) > 0) {
        for (name_idx in seq_along(sem_names)) {
          covariate_name <- sem_names[[name_idx]]
          effect_type <- .fsem_effect_type(
            res[[fit_key]]$mod$regression,
            latent_name,
            covariate_name,
            "response",
            "covariate"
          )
          denom <- if (identical(effect_type, "historical")) 200^2 else 200
          entry$coef.sem[[latent_name]][[covariate_name]] <- sum(
            (est_eval$coef.sem.std[[latent_name]][[name_idx]] -
               true_eval$coef.sem.std[[latent_name]][[name_idx]])^2
          ) / denom
        }
      }

      n_sem_cols <- ncol(est_eval$sigma.sem.eigvec.std[[sem_idx]])
      entry$sigma.sem.eigvec[[latent_name]] <- vapply(seq_len(n_sem_cols), function(col_idx) {
        mean((est_eval$sigma.sem.eigvec.std[[sem_idx]][, col_idx] -
                true_eval$sigma.sem.eigvec.std[[sem_idx]][, col_idx])^2)
      }, numeric(1))
      entry$sigma.sem.eigval[[latent_name]] <- vapply(seq_along(est_eval$sigma.sem.eigval.std[[sem_idx]]), function(col_idx) {
        (est_eval$sigma.sem.eigval.std[[sem_idx]][[col_idx]] -
           true_eval$sigma.sem.eigval.std[[sem_idx]][[col_idx]])^2
      }, numeric(1))
    }

    norm_list[[res_idx]] <- entry
  }

  mse <- vector("list", nrow(combo_df))
  for (combo_idx in seq_len(nrow(combo_df))) {
    N_target <- combo_df$N[[combo_idx]]
    M_target <- combo_df$M[[combo_idx]]
    matching <- Filter(function(entry) entry$N == N_target && entry$M == M_target, norm_list)
    if (length(matching) == 0) {
      next
    }

    combo_entry <- list(N = N_target, M = M_target, numb.runs = length(matching))

    for (indicator_name in nam.fac) {
      combo_entry[[indicator_name]] <- structure(list(), class = "data.frame", row.names = 1L)
      first_match <- matching[[1]]
      fac_names <- names(first_match$coef.fac[[indicator_name]])

      for (factor_name in fac_names) {
        combo_entry[[indicator_name]][[factor_name]] <- mean(vapply(matching, function(entry) {
          entry$coef.fac[[indicator_name]][[factor_name]]
        }, numeric(1)))
      }

      combo_entry[[indicator_name]]$intercept <- mean(vapply(matching, function(entry) {
        entry$intercept[[indicator_name]]
      }, numeric(1)))

      n_fac_cols <- length(first_match$sigma.fac.eigvec[[indicator_name]])
      for (col_idx in seq_len(n_fac_cols)) {
        combo_entry[[indicator_name]][[paste0("eig.fuc.", col_idx)]] <- mean(vapply(matching, function(entry) {
          entry$sigma.fac.eigvec[[indicator_name]][[col_idx]]
        }, numeric(1)))
        combo_entry[[indicator_name]][[paste0("eig.val.", col_idx)]] <- mean(vapply(matching, function(entry) {
          entry$sigma.fac.eigval[[indicator_name]][[col_idx]]
        }, numeric(1)))
      }

      combo_entry[[indicator_name]]$sig.err <- mean(vapply(matching, function(entry) {
        entry$sigma.error[[indicator_name]]
      }, numeric(1)))
    }

    for (latent_name in nam.sem) {
      combo_entry[[latent_name]] <- structure(list(), class = "data.frame", row.names = 1L)
      first_match <- matching[[1]]
      sem_names <- names(first_match$coef.sem[[latent_name]])

      for (covariate_name in sem_names) {
        combo_entry[[latent_name]][[covariate_name]] <- mean(vapply(matching, function(entry) {
          entry$coef.sem[[latent_name]][[covariate_name]]
        }, numeric(1)))
      }

      n_sem_cols <- length(first_match$sigma.sem.eigvec[[latent_name]])
      for (col_idx in seq_len(n_sem_cols)) {
        combo_entry[[latent_name]][[paste0("eig.fuc.", col_idx)]] <- mean(vapply(matching, function(entry) {
          entry$sigma.sem.eigvec[[latent_name]][[col_idx]]
        }, numeric(1)))
        combo_entry[[latent_name]][[paste0("eig.val.", col_idx)]] <- mean(vapply(matching, function(entry) {
          entry$sigma.sem.eigval[[latent_name]][[col_idx]]
        }, numeric(1)))
      }
    }

    mse[[combo_idx]] <- combo_entry
  }

  tabl.fac <- data.frame(N = combo_df$N, M = combo_df$M)
  for (combo_idx in seq_len(nrow(combo_df))) {
    for (fac_idx in seq_along(nam.fac)) {
      indicator_name <- nam.fac[[fac_idx]]
      indicator_mse <- mse[[combo_idx]][[indicator_name]]
      tabl.fac[[paste0("beta", fac_idx)]][combo_idx] <- indicator_mse$intercept

      fac_names <- names(estimation$result$params.estimated.eval$coef.fac.std[[fac_idx]])
      for (name_idx in seq_along(fac_names)) {
        tabl.fac[[paste0("lambda", fac_idx, name_idx)]][combo_idx] <- indicator_mse[[fac_names[[name_idx]]]]
      }

      n_fac_cols <- ncol(estimation$result$params.estimated.eval$sigma.fac.eigvec.std[[fac_idx]])
      for (col_idx in seq_len(n_fac_cols)) {
        tabl.fac[[paste0("phi", fac_idx, col_idx)]][combo_idx] <- indicator_mse[[paste0("eig.fuc.", col_idx)]]
        tabl.fac[[paste0("nu", fac_idx, col_idx)]][combo_idx] <- indicator_mse[[paste0("eig.val.", col_idx)]]
      }

      tabl.fac[[paste0("sigma", fac_idx)]][combo_idx] <- indicator_mse$sig.err
    }
  }

  tabl.sem <- data.frame(N = combo_df$N, M = combo_df$M)
  for (combo_idx in seq_len(nrow(combo_df))) {
    for (sem_idx in seq_along(nam.sem)) {
      latent_name <- nam.sem[[sem_idx]]
      latent_mse <- mse[[combo_idx]][[latent_name]]
      sem_names <- names(estimation$result$params.estimated.eval$coef.sem.std[[latent_name]])

      for (name_idx in seq_along(sem_names)) {
        tabl.sem[[paste0("gamma", sem_idx, name_idx)]][combo_idx] <- latent_mse[[sem_names[[name_idx]]]]
      }

      n_sem_cols <- ncol(estimation$result$params.estimated.eval$sigma.sem.eigvec.std[[sem_idx]])
      for (col_idx in seq_len(n_sem_cols)) {
        tabl.sem[[paste0("psi", sem_idx, col_idx)]][combo_idx] <- latent_mse[[paste0("eig.fuc.", col_idx)]]
        tabl.sem[[paste0("mu", sem_idx, col_idx)]][combo_idx] <- latent_mse[[paste0("eig.val.", col_idx)]]
      }
    }
  }

  list(mse = mse, tabl.fac = tabl.fac, tabl.sem = tabl.sem)
}

if (!exists("results", inherits = TRUE)) {
  stop("results object not found")
}

if (length(results) == 0) {
  stop("results is empty")
}

if (!is.null(results[[1]]$estimation)) {
  mse_report <- .fsem_compute_single_mse(results, "model.fit", "estimation")
  mse <- mse_report$mse
  tabl.fac <- mse_report$tabl.fac
  tabl.sem <- mse_report$tabl.sem
} else if (!is.null(results[[1]]$estimation.w) && !is.null(results[[1]]$estimation.m)) {
  mse_report.w <- .fsem_compute_single_mse(results, "model.fit.w", "estimation.w")
  mse_report.m <- .fsem_compute_single_mse(results, "model.fit.m", "estimation.m")
  mse.w <- mse_report.w$mse
  mse.m <- mse_report.m$mse
  tabl.fac.w <- mse_report.w$tabl.fac
  tabl.sem.w <- mse_report.w$tabl.sem
  tabl.fac.m <- mse_report.m$tabl.fac
  tabl.sem.m <- mse_report.m$tabl.sem
} else {
  stop("results structure is not supported by Core/MSE_report.R")
}