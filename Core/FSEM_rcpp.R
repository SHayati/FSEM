######################################################################
## FSEM_rcpp.R
##
## em.estimation.rcpp: a further-accelerated re-implementation of
## em.estimation() that pushes the inner CV + sig-iter loops (which are the
## dominant cost after the closed-form penalty optimisation) entirely into
## Rcpp via Core/fsem_kernels_rcpp.cpp.
##
## The ORIGINAL implementation (em.estimation in Core/FSEM.R) is left
## untouched. This file ALSO does not touch em.estimation.optimised; it
## reuses its penalty-precompute helpers (.fsem_precompute_penalties etc.).
##
## Statistical algorithm is identical to em.estimation / em.estimation.optimised.
## Verified to machine precision against the original at n.em = 1.
##
## Requires: source Core/FSEM.R first, then Core/FSEM_optimised.R, then this file.
######################################################################

.fsem_rcpp_env  <- new.env(parent = globalenv())
.fsem_have_rcpp_full <- NA

.fsem_load_rcpp_kernels <- function() {
  if (!is.na(.fsem_have_rcpp_full)) return(invisible(.fsem_have_rcpp_full))
  ok <- FALSE
  if (requireNamespace("Rcpp", quietly = TRUE)) {
    candidates <- c(
      file.path("Core", "fsem_kernels_rcpp.cpp"),
      "fsem_kernels_rcpp.cpp",
      file.path(getwd(), "Core", "fsem_kernels_rcpp.cpp")
    )
    kp <- candidates[file.exists(candidates)][1]
    if (!is.na(kp)) {
      ok <- tryCatch({
        old_libs <- Sys.getenv("PKG_LIBS", unset = NA)
        Sys.setenv(PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)")
        on.exit({
          if (is.na(old_libs)) Sys.unsetenv("PKG_LIBS")
          else Sys.setenv(PKG_LIBS = old_libs)
        }, add = TRUE)
        Rcpp::sourceCpp(kp, env = .fsem_rcpp_env, rebuild = TRUE)
        TRUE
      }, error = function(e) {
        message("em.estimation.rcpp: Rcpp kernels unavailable (",
                conditionMessage(e), "); falling back to em.estimation.optimised.")
        FALSE
      })
    }
  }
  assign(".fsem_have_rcpp_full", ok, envir = topenv())
  invisible(ok)
}

## Build per-sample design matrices f[[m]][[k]] and f2[[m]][[k]] for a factor j.
## Identical algebra to em.estimation.optimised().
.fsem_build_fac_design <- function(model, j, samples, mont, no.lat, n.b,
                                   eta.l, eval.t, omegaMatrix,
                                   design.regular, soll_global) {
  f  <- vector("list", mont)
  f2 <- vector("list", mont)
  for (m in 1:mont) {
    f[[m]]  <- vector("list", samples); names(f[[m]])  <- as.character(1:samples)
    f2[[m]] <- vector("list", samples); names(f2[[m]]) <- as.character(1:samples)
    for (k in 1:samples) {
      soll <- if (design.regular) soll_global
              else ginv(eval.t$eval[[k]][[j]] %*% t(eval.t$eval[[k]][[j]])) %*% eval.t$eval[[k]][[j]]
      Fk <- diag(n.b)
      f2k <- 0
      for (i2 in 1:no.lat) {
        et <- kronecker(diag(eval.t$time.no[[k]][[j]]), t(eta.l[[m]][[k]][[i2]]))
        fac <- which(model$mod$factorModel[[j]]$factor == model$var$latents[i2])
        if (length(fac) > 0) {
          eff <- model$mod$factorModel[[j]]$effect[fac]
          ef <- NULL
          if (eff == "concurrent") { ef <- soll %*% (et %*% omegaMatrix$omega1[[k]][[j]]); f2k <- 0 }
          else if (eff == "historical") { ef <- soll %*% (et %*% omegaMatrix$omega2[[k]][[j]]); f2k <- 0 }
          else if (eff == "fixed_concurrent") { ef <- NULL; f2k <- eta.l[[m]][[k]][[i2]] }
          else if (eff == "fixed_historical") { ef <- NULL; f2k <- soll %*% omegaMatrix$delta[[k]][[j]] %*% eta.l[[m]][[k]][[i2]] }
          if (!is.null(ef)) Fk <- cbind(Fk, ef)
        }
      }
      f[[m]][[as.character(k)]]  <- Fk
      f2[[m]][[as.character(k)]] <- f2k
    }
  }
  list(f = f, f2 = f2)
}

## Stack per-sample matrices/vectors for a chosen sample subset.
## Returns Fs (list of M stacked matrices) and rs (list of M residual vectors).
.fsem_stack_fac <- function(samp, f, f2, gen_z, j, mont, n.b) {
  Fs <- vector("list", mont)
  rs <- vector("list", mont)
  for (m in 1:mont) {
    Fs[[m]] <- do.call(rbind, lapply(samp, function(i1) f[[m]][[as.character(i1)]]))
    zz <- unlist(lapply(samp, function(i1) gen_z[[as.character(i1)]][[j]][[m]]))
    ff2 <- unlist(lapply(samp, function(i1) {
      x <- f2[[m]][[as.character(i1)]]
      if (length(x) == 1L) rep(x, n.b) else as.vector(x)
    }))
    if (length(ff2) != length(zz)) {
      ff2 <- rep_len(ff2, length(zz))   # match R recycling behaviour
    }
    rs[[m]] <- as.vector(zz) - ff2
  }
  list(Fs = Fs, rs = rs)
}

## Stack per-sample matrices/vectors for SEM. ett[[m]] is response.
## g is the per-sample design (list of M of list of samples).
.fsem_stack_sem <- function(samp, g, eta.l, lat.count, mont) {
  Fs <- vector("list", mont)
  rs <- vector("list", mont)
  for (m in 1:mont) {
    Fs[[m]] <- do.call(rbind, lapply(samp, function(i1) g[[m]][[i1]]))
    rs[[m]] <- unlist(lapply(samp, function(i1) eta.l[[m]][[i1]][[lat.count]]))
  }
  list(Fs = Fs, rs = rs)
}

######################################################################
## em.estimation.rcpp : Rcpp-fused version of em.estimation.optimised
######################################################################
em.estimation.rcpp <- function(model, data, x.data = NULL,
                                initial.parameter, n.b, n.em = 100, n.monte = 100,
                                s.p = c("min", "sd"), range.min = NULL, range.max = NULL,
                                design = c("regular", "irregular", "regular.truncated", "regular.missing"),
                                plot.progress = FALSE,
                                save.every = NULL) {

  s.p = match.arg(s.p)
  if (is.null(range.min)) range.min = 0
  if (is.null(range.max)) range.max = 1
  if (is.null(s.p)) s.p = "min"

  ok <- .fsem_load_rcpp_kernels()
  if (!isTRUE(ok)) {
    message("em.estimation.rcpp: kernels not loaded; delegating to em.estimation.optimised().")
    return(em.estimation.optimised(model, data, x.data, initial.parameter, n.b,
                                   n.em, n.monte, s.p, range.min, range.max,
                                   design, plot.progress, save.every))
  }
  fit_inner <- .fsem_rcpp_env$fsem_fit_inner
  cv_fit    <- .fsem_rcpp_env$fsem_cv_fit

  omegaMatrix <- omega.matrix(n.b, range.min, range.max, data)
  result <- list()
  no.em <- n.em
  no.it <- 10
  parameterr <- initial.parameter
  diagMatrix  <- diag.matrix(model, n.b)
  diagMatrix2 <- diag.matrix2(model, n.b)
  samples <- length(unique(data$.id))
  weightMatrix <- weight.matrix(model, n.b, sample = samples)
  no.fac <- length(model$mod$factorModel)
  no.lat <- length(model$var$latents)
  basis <- create.bspline.basis(nbasis = n.b, rangeval = c(range.min, range.max))
  pen <- getbasispenalty(basis)
  ome <- inprod(basis, basis)
  eval.t <- timeEval(n.b, range.min, range.max, data)
  n.t <- max(c(sapply(1:samples, function(i) { eval.t$time.no[[i]] })))
  evall2 <- t(eval.basis(seq(0, 1, length.out = n.t), basis))
  mont <- n.monte
  design = match.arg(design)
  design.regular <- (design == "regular")

  ## Precompute constant penalty matrices (closed form, from FSEM_optimised.R)
  penalties <- .fsem_precompute_penalties(model, basis, n.b, x.data, range.min, range.max)
  pen.fac.list <- penalties$fac
  pen.sem.list <- penalties$sem

  pb <- txtProgressBar(0, no.em, style = 3)
  pbk <- 0

  res.list <- list()
  cond <- list()
  res.k <- 0

  likeli.fac <- c(); penal.fac <- c(); likeli.sem <- c(); penal.sem <- c()
  toler_mean <- c(); toler_mean_sig <- c()

  for (i in 1:no.em) {

    eta.distt <- eta.distribution(model, param = parameterr, n.b, range.min, range.max, data, x.data, omegaMatrix = omegaMatrix)
    paramz1 <- params.z1(model, param = parameterr, n.b, range.min, range.max, eta.dist = eta.distt, data, omegaMatrix = omegaMatrix, design = design)
    time.start <- Sys.time()

    ################ Monte-Carlo z draws (same as original) ##################
    generators2.z <- list()
    generators2.w <- list()
    for (i1 in 1:samples) {
      samp <- i1
      gg <- list()
      mu.gg <- list()
      generators2.z[[paste0(samp)]] <- list()
      generators2.w[[paste0(samp)]] <- list()
      w.tot <- data$.value[data$.id == samp]
      paramz <- params.z(model, param = paramz1, n.b = n.b,
                         range.min = range.min, range.max = range.max,
                         sample.no = i1, w = w.tot, data)
      for (j1 in 1:no.fac) {
        generators2.z[[paste0(samp)]][[j1]] <- list()
        generators2.w[[paste0(samp)]][[j1]] <- list()
        generators3.w <- data$.value[data$.ind == j1 & data$.id == samp]
        for (m in 1:mont) {
          if (j1 == 1) {
            generators.z <- generator.z(model, param = paramz, indicator.no = j1)
            gg[[j1]] <- generators.z$zw.tot
            mu.gg[[j1]] <- as.vector(paramz$mu.zw.j[[j1]])
          } else {
            generators.z <- generator.z(model, param = paramz, indicator.no = j1,
                                        zw = gg[[j1 - 1]], mu.zw = mu.gg[[j1 - 1]])
            gg[[j1]] <- c(as.vector(gg[[j1 - 1]]), as.vector(generators.z$zw.tot))
            mu.gg[[j1]] <- c(mu.gg[[j1 - 1]], as.vector(paramz$mu.zw.j[[j1]]))
          }
          generators2.z[[paste0(samp)]][[j1]][[m]] <- generators.z$zw.tot
          generators2.w[[paste0(samp)]][[j1]][[m]] <- generators3.w
        }
      }
    }

    beta.s <- list()
    for (m in 1:mont) {
      beta.s[[m]] <- list()
      for (j1 in 1:no.lat) {
        beta.s[[m]][[j1]] <- list()
        for (i1 in 1:samples) beta.s[[m]][[j1]][[i1]] <- Inf
      }
    }

    paramzz <- paramz1
    xi <- list()
    for (i1 in 1:samples) {
      samp <- i1
      xi[[paste0(samp)]] <- list()
      zz <- c(sapply(1:no.fac, function(j1) { as.vector(rowMeans(sapply(1:mont, function(ii) { generators2.z[[paste0(samp)]][[j1]][[ii]] }))) }))
      for (m in 1:mont) {
        ## NOTE: kept arithmetically identical to original to preserve eigen
        ## sign / RNG realisations in mvrnorm().
        mu.eta.z <- t(t(paramzz[[i1]]$params$mu.eta)) + paramzz[[i1]]$params$sigma.eta %*%
          t(paramzz[[i1]]$params$a.matrix) %*% solve(paramzz[[i1]]$params$sigma.z) %*%
          t(t(zz - paramzz[[i1]]$params$mean.z))
        sigma.eta.z <- paramzz[[i1]]$params$sigma.eta - paramzz[[i1]]$params$sigma.eta %*%
          t(paramzz[[i1]]$params$a.matrix) %*% solve(paramzz[[i1]]$params$sigma.z) %*%
          paramzz[[i1]]$params$a.matrix %*% paramzz[[i1]]$params$sigma.eta
        sigma.eta.z <- 0.5 * (sigma.eta.z + t(sigma.eta.z))
        loop.count <- 0
        repeat {
          eta.zz <- as.vector(mvrnorm(1, mu.eta.z, sigma.eta.z))
          a <- c(); count2 <- 0
          for (i2 in 1:no.lat) {
            etaa <- eta.zz[(count2 + 1):(count2 + n.b)]
            count2 <- count2 + n.b
            a[i2] <- t(etaa) %*% pen %*% etaa <= beta.s[[m]][[i2]][[samp]]
          }
          loop.count <- loop.count + 1
          if (all(a)) break
          if (loop.count > 20) break
        }
        xi[[paste0(samp)]][[m]] <- eta.zz
      }
    }

    eta.l <- list()
    for (m in 1:mont) {
      eta.l[[m]] <- list()
      for (k in 1:samples) {
        eta.l[[m]][[k]] <- list()
        count2 <- 0
        for (i2 in 1:no.lat) {
          eta.l[[m]][[k]][[i2]] <- xi[[paste0(k)]][[m]][(count2 + 1):(count2 + n.b)]
          count2 <- count2 + n.b
        }
      }
    }
    for (m in 1:mont) {
      for (i2 in 1:no.lat) {
        meanv <- rowMeans(sapply(1:samples, function(k) { eta.l[[m]][[k]][[i2]] }))
        for (k in 1:samples) eta.l[[m]][[k]][[i2]] <- eta.l[[m]][[k]][[i2]] - meanv
      }
    }

    if (design.regular == TRUE) {
      soll <- ginv(eval.t$eval[[1]][[1]] %*% t(eval.t$eval[[1]][[1]])) %*% eval.t$eval[[1]][[1]]
    } else {
      soll <- NULL
    }

    ################ factor parameter estimations ############################
    delta.fac <- c()
    hist.fac <- list()
    for (j in 1:no.fac) {
      aa <- c()
      for (k in 1:dim(diagMatrix2[[paste0("f", j)]])[1]) {
        aa[k] <- if (all(diagMatrix2[[paste0("f", j)]][k, ] == 0)) 0 else 1
      }
      pen.fac <- pen.fac.list[[j]]

      des <- .fsem_build_fac_design(model, j, samples, mont, no.lat, n.b,
                                    eta.l, eval.t, omegaMatrix,
                                    design.regular, soll)
      f  <- des$f
      f2 <- des$f2

      if (i > 0) {
        delt_sta <- c(0.001, 0.01, 0.05, 0.1) * log(samples)
        delt <- delt_sta
        samp.train <- sort(sample(1:samples, round(0.7 * samples)))
        samp.test  <- sort((1:samples)[-samp.train])

        train <- .fsem_stack_fac(samp.train, f, f2, generators2.z, j, mont, n.b)
        test  <- .fsem_stack_fac(samp.test , f, f2, generators2.z, j, mont, n.b)

        sig <- parameterr$sigma.fac[[j]]
        cv <- cv_fit(train$Fs, train$rs, test$Fs, test$rs,
                     sig, pen.fac, as.numeric(delt), as.integer(no.it), as.integer(n.b))
        per       <- as.numeric(cv$per)
        lambda2.1 <- matrix(as.numeric(cv$lambda), ncol = 1)
        ## sig from CV (last delta) intentionally not used; original final fit
        ## resets to parameterr$sigma.fac[[j]] below.

        delt.min <- delt[which(per == min(per))]
        ex1 <- vapply(seq_len(mont), function(m) {
          A <- cv$A_test[[m]]; b <- cv$b_test[[m]]
          as.numeric((-t(lambda2.1) %*% b) +
                     t(lambda2.1) %*% (A + delt.min * pen.fac) %*% lambda2.1)
        }, numeric(1))
        bound <- mean(ex1) + sd(ex1)
        per2 <- vapply(seq_along(delt), function(l) {
          v <- 0
          for (m in 1:mont) {
            A <- cv$A_test[[m]]; b <- cv$b_test[[m]]
            v <- v + 1/mont * as.numeric((-t(lambda2.1) %*% b) +
                                          t(lambda2.1) %*% (A + delt[l] * pen.fac) %*% lambda2.1)
          }
          as.numeric(v <= bound)
        }, numeric(1))
        d.t <- delt[as.logical(per2)]
        deltaa <- if (s.p == "min") delt.min else max(d.t)
      } else {
        deltaa <- ifelse(n.b > 2, 0.01, 0)
      }
      delta.fac[j] <- deltaa

      ## Final fit on ALL samples with chosen delta
      full <- .fsem_stack_fac(1:samples, f, f2, generators2.z, j, mont, n.b)
      sig <- parameterr$sigma.fac[[j]]
      fit <- fit_inner(full$Fs, full$rs, sig, pen.fac,
                       as.numeric(deltaa), as.integer(no.it), as.integer(n.b))
      lambda2.1 <- matrix(as.numeric(fit$lambda), ncol = 1)
      sig <- fit$sig

      lambda2 <- c(); count <- 1
      for (k in 1:length(aa)) {
        if (aa[k] == 1) { lambda2[k] <- lambda2.1[count]; count <- count + 1 } else lambda2[k] <- 0
      }
      sig.err <- 0
      for (i1 in 1:samples) {
        zsum <- 0
        for (m in 1:mont) {
          sigg <- generators2.w[[paste0(i1)]][[j]][[1]] - t(eval.t$eval[[i1]][[j]]) %*% generators2.z[[paste0(i1)]][[j]][[m]]
          zsum <- zsum + 1 / mont * (t(sigg) %*% sigg)
        }
        sig.err <- sig.err + zsum
      }
      sig.err1 <- as.vector(sig.err / sum(sapply(1:samples, function(i1) eval.t$time.no[[i1]][[j]])))
      parameterr$lambda.param[[j]] <- as.vector(lambda2)
      parameterr$sigma.fac[[j]] <- sig
      parameterr$sigma.error[[j]] <- sig.err1
      hist.fac[[j]] <- list(generators2.z = generators2.z, generators2.w = generators2.w,
                            eta.l = eta.l, f = f, f2 = f2, deltaa = deltaa, pen.fac = pen.fac)
    }

    ################ SEM parameter estimations ###############################
    no.reg <- length(model$mod$regression)
    lat2 <- model$var$latents
    obs <- model$var$observed
    no.lat2 <- length(lat2)
    no.obs2 <- length(obs)
    soll <- ginv(evall2 %*% t(evall2)) %*% evall2
    hist.sem <- list()
    for (j in 1:no.reg) {
      co <- model$mod$regression[[j]]$covariate
      no.cov <- length(co)
      lat.count <- which(lat2 == model$mod$regression[[j]]$response)
      aa <- c()
      if (length(obs) > 0) {
        ddm <- Matrix::bdiag(diagMatrix2[[paste0("r", j, ".eta")]], diagMatrix2[[paste0("r", j, ".x")]])
        dd <- dim(ddm)
        for (kk in 1:dd[1]) aa[kk] <- if (all(ddm[kk, ] == 0)) 0 else 1
      } else {
        dd <- dim(diagMatrix2[[paste0("r", j, ".eta")]])
        for (kk in 1:dd[1]) aa[kk] <- if (all(diagMatrix2[[paste0("r", j, ".eta")]][kk, ] == 0)) 0 else 1
      }
      if (length(which(aa == 1)) > 0) {
        pen.sem <- pen.sem.list[[j]]
      }
      g.eta <- NULL; g.x <- NULL
      et.ex_present <- FALSE
      x.ex_present  <- FALSE
      if (no.cov > 0) {
        if (any(co %in% model$var$latents)) {
          g.eta <- list()
          for (m in 1:mont) {
            g.eta[[m]] <- list()
            for (i1 in 1:samples) {
              g.eta.tmp <- NULL
              for (jj in 1:no.lat2) {
                fac <- which(model$mod$regression[[j]]$covariate == lat2[jj])
                if (length(fac) == 0) {
                  gg.eta <- NULL
                } else {
                  if (model$mod$regression[[j]]$effect[fac] == "concurrent") {
                    gg.eta <- soll %*% kronecker(diag(n.t), t(eta.l[[m]][[i1]][[jj]])) %*% omegaMatrix$omega1.sem
                  }
                  if (model$mod$regression[[j]]$effect[fac] == "historical") {
                    gg.eta <- soll %*% kronecker(diag(n.t), t(eta.l[[m]][[i1]][[jj]])) %*% omegaMatrix$omega2.sem
                  }
                }
                g.eta.tmp <- cbind(g.eta.tmp, gg.eta)
              }
              g.eta[[m]][[i1]] <- g.eta.tmp
            }
          }
        } else {
          g.eta <- lapply(1:mont, function(m) { lapply(1:samples, function(i1) (g.eta[[m]][[i1]] = 0)) })
          et.ex_present <- TRUE; et.ex <- "NULL"
        }
        if (length(obs) > 0) {
          if (any(co %in% model$var$observed)) {
            g.x <- list(); xx <- list()
            for (jj in 1:no.obs2) {
              fac <- which(model$mod$regression[[j]]$covariate == obs[jj])
              if (model$mod$regression[[j]]$effect[fac] %in% c("concurrent", "historical")) {
                sm <- matrix(NA, n.b, samples)
                for (tt in 1:samples) {
                  sm[, tt] <- Data2fd(argvals = seq(range.min, range.max, length.out = length(x.data[[obs[jj]]][[tt]])), y = x.data[[obs[jj]]][[tt]], basisobj = basis, lambda = 0.5)$coefs
                }
                meann <- rowMeans(sm); xx[[jj]] <- sm - meann
              }
              if (model$mod$regression[[j]]$effect[fac] == "smooth") {
                basis3 <- create.bspline.basis(rangeval = range(x.data[[model$var$observed[jj]]]), nbasis = n.b)
                sm <- t(eval.basis(x.data[[model$var$observed[jj]]], basis3))
                meann <- rowMeans(sm); xx[[jj]] <- sm - meann
              }
              if (model$mod$regression[[j]]$effect[fac] == "linear") {
                if (length(unique(x.data[[model$var$observed[jj]]])) == 2) {
                  xx[[jj]] <- x.data[[model$var$observed[jj]]]
                } else {
                  meann <- mean(x.data[[model$var$observed[jj]]]); xx[[jj]] <- x.data[[model$var$observed[jj]]] - meann
                }
              }
            }
            for (i1 in 1:samples) {
              g.x.tmp <- NULL
              for (jj in 1:no.obs2) {
                fac <- which(model$mod$regression[[j]]$covariate == obs[jj])
                if (length(fac) == 0) {
                  gg.x <- NULL
                } else {
                  if (model$mod$regression[[j]]$effect[fac] == "concurrent") {
                    gg.x <- soll %*% kronecker(diag(n.t), t(xx[[jj]][, i1])) %*% omegaMatrix$omega1.sem
                  }
                  if (model$mod$regression[[j]]$effect[fac] == "historical") {
                    gg.x <- soll %*% kronecker(diag(n.t), t(xx[[jj]][, i1])) %*% omegaMatrix$omega2.sem
                  }
                  if (model$mod$regression[[j]]$effect[fac] == "smooth") {
                    gg.x <- soll %*% t(kronecker(evall2, xx[[jj]][, i1]))
                  }
                  if (model$mod$regression[[j]]$effect[fac] == "linear") {
                    gg.x <- as.vector(xx[[jj]][i1]) * diag(n.b)
                  }
                }
                g.x.tmp <- cbind(g.x.tmp, gg.x)
                g.x[[i1]] <- g.x.tmp
              }
            }
          } else {
            g.x <- lapply(1:samples, function(i1) (g.x[[i1]] = 0))
            x.ex_present <- TRUE; x.ex <- "NULL"
          }
        }
      }

      sig <- parameterr$sigma.sem[[lat.count]]
      if (!is.null(g.eta) | !is.null(g.x)) {
        if (i > 0) {
          delt_sta <- c(0.001, 0.01, 0.05, 0.1) * log(samples)
          delt <- delt_sta
          samp.train <- sort(sample(1:samples, round(0.7 * samples)))
          samp.test  <- sort((1:samples)[-samp.train])

          ## Build per-sample g[[m]][[i1]] just like the original
          g <- vector("list", mont)
          for (m in 1:mont) {
            g[[m]] <- vector("list", samples)
            for (i1 in 1:samples) {
              gi <- cbind(g.eta[[m]][[i1]], g.x[[i1]])
              if (et.ex_present) gi <- gi[, -1, drop = FALSE]
              if (x.ex_present)  gi <- gi[, -ncol(gi), drop = FALSE]
              g[[m]][[i1]] <- gi
            }
          }

          train <- .fsem_stack_sem(samp.train, g, eta.l, lat.count, mont)
          test  <- .fsem_stack_sem(samp.test , g, eta.l, lat.count, mont)

          cv <- cv_fit(train$Fs, train$rs, test$Fs, test$rs,
                       sig, pen.sem, as.numeric(delt),
                       as.integer(no.it), as.integer(n.b))
          per       <- as.numeric(cv$per)
          gamma2.1 <- matrix(as.numeric(cv$lambda), ncol = 1)

          delt.min <- delt[which(per == min(per))]
          ex1 <- vapply(seq_len(mont), function(m) {
            A <- cv$A_test[[m]]; b <- cv$b_test[[m]]
            as.numeric((-t(gamma2.1) %*% b) +
                       t(gamma2.1) %*% (A + delt.min * pen.sem) %*% gamma2.1)
          }, numeric(1))
          bound <- mean(ex1) + sd(ex1)
          per2 <- vapply(seq_along(delt), function(l) {
            v <- 0
            for (m in 1:mont) {
              A <- cv$A_test[[m]]; b <- cv$b_test[[m]]
              v <- v + 1/mont * as.numeric((-t(gamma2.1) %*% b) +
                                            t(gamma2.1) %*% (A + delt[l] * pen.sem) %*% gamma2.1)
            }
            as.numeric(v <= bound)
          }, numeric(1))
          d.t <- delt[as.logical(per2)]
          deltaa <- if (s.p == "min") delt.min else max(d.t)
        } else {
          deltaa <- 0
        }
        delta.sem <- deltaa

        ## Final fit on ALL samples with chosen delta
        g <- vector("list", mont)
        for (m in 1:mont) {
          g[[m]] <- vector("list", samples)
          for (i1 in 1:samples) {
            gi <- cbind(g.eta[[m]][[i1]], g.x[[i1]])
            if (et.ex_present) gi <- gi[, -1, drop = FALSE]
            if (x.ex_present)  gi <- gi[, -ncol(gi), drop = FALSE]
            g[[m]][[i1]] <- gi
          }
        }
        full <- .fsem_stack_sem(1:samples, g, eta.l, lat.count, mont)
        sig <- parameterr$sigma.sem[[lat.count]]
        fit <- fit_inner(full$Fs, full$rs, sig, pen.sem,
                         as.numeric(deltaa), as.integer(no.it), as.integer(n.b))
        gamma2.1 <- matrix(as.numeric(fit$lambda), ncol = 1)
        sig <- fit$sig
      } else {
        ## No covariates: sigma = mean over samples and m of eta.l_outer
        sigma1 <- 0
        for (i1 in 1:samples) {
          ss2 <- 0
          for (m in 1:mont) {
            ss <- eta.l[[m]][[i1]][[lat.count]]
            ss2 <- ss2 + 1 / mont * ss %*% t(ss)
          }
          sigma1 <- sigma1 + 1 / samples * ss2
        }
        sig <- sigma1
      }

      gamma2 <- c(); count <- 1
      for (k in 1:length(aa)) {
        if (aa[k] == 1) { gamma2[k] <- gamma2.1[count]; count <- count + 1 } else gamma2[k] <- 0
      }
      parameterr$gamma.param[[lat.count]] <- as.vector(gamma2)
      parameterr$sigma.sem[[lat.count]] <- sig
      if (length(which(aa == 1)) > 0) {
        hist.sem[[lat.count]] <- list(g.eta = g.eta, g.x = g.x, g = if (exists("g")) g else NULL,
                                      deltaa = deltaa, pen.sem = pen.sem)
      } else {
        hist.sem[[lat.count]] <- list(g.eta = g.eta, g.x = g.x, deltaa = deltaa)
      }
    }

    pbk <- pbk + 1
    setTxtProgressBar(pb, pbk)
    time.end <- Sys.time()
    parameters <- ex.estimations(model, estimate = parameterr, n.b = n.b, range.min = range.min, range.max = range.max, data = data, x.data = x.data)
    parameterr <- parameters$totparam
    params.eval <- parameter.estimated.evaluted(model, param = parameters, n.b = n.b, n.monte = n.monte, range.min = range.min, range.max = range.max, data = data, x.data = x.data)

    ### convergence ###
    if (i > 1) {
      count <- 0; toll <- c()
      for (j in 1:no.fac) {
        ind <- model$mod$factorModel[[j]]$indicator
        for (jj in 1:length(model$mod$factorModel[[j]]$factor)) {
          fac <- model$mod$factorModel[[j]]$factor[[jj]]
          par <- params.eval$coef.fac.std[[paste0(ind)]][[paste0(fac)]]
          par.p <- param_prev$coef.fac.std[[paste0(ind)]][[paste0(fac)]]
          rim <- sqrt(mean((par - par.p)^2)); count <- count + 1; toll[count] <- rim
        }
        par <- params.eval$intercept.std[[paste0(ind)]]
        par.p <- param_prev$intercept.std[[paste0(ind)]]
        rim <- sqrt(mean((par - par.p)^2)); count <- count + 1; toll[count] <- rim
      }
      for (j in 1:length(model$mod$regression)) {
        res <- model$mod$regression[[j]]$response
        if (length(model$mod$regression[[j]]$covariate) > 0) {
          for (jj in 1:length(model$mod$regression[[j]]$covariate)) {
            cov <- model$mod$regression[[j]]$covariate[[jj]]
            par <- params.eval$coef.sem.std[[paste0(res)]][[paste0(cov)]]
            par.p <- param_prev$coef.sem.std[[paste0(res)]][[paste0(cov)]]
            rim <- sqrt(mean((par - par.p)^2)); count <- count + 1; toll[count] <- rim
          }
        }
      }
      toler_mean[i] <- max(toll)
      toll_sig <- c()
      for (j in 1:length(model$mod$factorModel)) {
        toll_sig[j] <- abs(parameterr$sigma.error[[j]] - sig_prev[[j]])
      }
      toler_mean_sig[i] <- max(toll_sig)
    }

    param_prev <- params.eval
    sig_prev <- parameterr$sigma.error
    res.k = res.k + 1
    res.list[[res.k]] = list(parameters = parameters, parameters.eval = params.eval, timing = c(time.start = time.start, time.end = time.end))

    if (!is.null(save.every)) {
      if (!is.null(save.every$steps) && (i %% save.every$steps) == 0) {
        file_name = gsub("[step]", i, tolower(save.every$file_name_format), fixed = TRUE)
        file_path = file.path(dir_path, file_name)
        result2 = result
        result2$params.estimated <- parameters
        result2$params.estimated.eval <- params.eval
        res_to_save <- list(result = result2,
                            tolerance = list(tol = toler_mean, tol_sig = toler_mean_sig),
                            steps = i,
                            hist = list(hist.fac = hist.fac, hist.sem = hist.sem, eta = eta.l),
                            attached_info = save.every$attached_info)
        saveRDS(res_to_save, file_path)
      }
    }
  }
  result$params.estimated <- parameters
  result$params.estimated.eval <- params.eval
  final <- list(result = result,
                tolerance = list(tol = toler_mean, tol_sig = toler_mean_sig), steps = i,
                hist = list(hist.fac = hist.fac, hist.sem = hist.sem, eta = eta.l))
  final
}
