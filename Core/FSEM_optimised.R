######################################################################
## FSEM_optimised.R
##
## Optimised re-implementation of em.estimation() from Core/FSEM.R.
##
## The ORIGINAL implementation in Core/FSEM.R is left untouched. This file
## defines em.estimation.optimised() with the SAME signature and the SAME
## statistical algorithm, but with the following performance optimisations:
##
##  1. Penalty matrices (the historical P+Ps+Pt and the smooth penalty), which
##     the original rebuilds inside the EM loop with a 100x100 double loop of
##     rank-1 updates on (n.b^2 x n.b^2) matrices, are computed ONCE before the
##     EM loop using an exact closed form derived from the Kronecker identity
##         kron(a,b) kron(a,b)^T = kron(a a^T, b b^T)
##     and the separability of the double sum:
##         P  = kron(E1 E1^T, E1 E1^T) / T^2
##         Ps = kron(E2 E2^T, E  E^T ) / T^2
##         Pt = kron(E  E^T , E2 E2^T) / T^2
##         Psmooth = kron(E2 E2^T, sm sm^T) / T
##     This is numerically identical to the original (verified to ~1e-12).
##
##  2. Every solve(kronecker(diag(n), sig)) is replaced by block operations:
##     since kron(diag(n), sig) is block-diagonal with n identical blocks,
##     its inverse is kron(diag(n), solve(sig)), and the quadratic / bilinear
##     forms t(F) %*% kron(I_n, Si) %*% F  and  t(F) %*% kron(I_n, Si) %*% z
##     are computed block-wise (O(n * nb^2 * p)) instead of inverting an
##     (n*nb) x (n*nb) matrix (O((n*nb)^3)).
##
##  3. Optional Rcpp kernels (Core/fsem_kernels.cpp) accelerate the block
##     forms further. If the kernels cannot be compiled/loaded, an equivalent
##     pure-R block implementation is used automatically.
##
## Requirements: source Core/FSEM.R first (this file reuses its helper
## functions: omega.matrix, diag.matrix, diag.matrix2, weight.matrix,
## timeEval, eta.distribution, params.z1, params.z, generator.z,
## ex.estimations, parameter.estimated.evaluted).
######################################################################

## ---- optional Rcpp kernels -----------------------------------------------
.fsem_kernel_env <- new.env(parent = globalenv())
.fsem_have_rcpp  <- NA  # NA = not yet attempted

.fsem_load_kernels <- function() {
  if (!is.na(.fsem_have_rcpp)) return(invisible(.fsem_have_rcpp))
  ok <- FALSE
  if (requireNamespace("Rcpp", quietly = TRUE)) {
    candidates <- c(
      file.path("Core", "fsem_kernels.cpp"),
      "fsem_kernels.cpp",
      file.path(getwd(), "Core", "fsem_kernels.cpp")
    )
    kp <- candidates[file.exists(candidates)][1]
    if (!is.na(kp)) {
      ok <- tryCatch({
        Rcpp::sourceCpp(kp, env = .fsem_kernel_env)
        TRUE
      }, error = function(e) {
        message("em.estimation.optimised: Rcpp kernels unavailable (",
                conditionMessage(e), "); using pure-R fallback.")
        FALSE
      })
    }
  }
  assign(".fsem_have_rcpp", ok, envir = topenv())
  invisible(ok)
}

## t(F) %*% kron(I_n, Si) %*% F   (F: (n*nb) x p, blocks of nb rows)
.fsem_xtSx <- function(F, Si, nb) {
  if (isTRUE(.fsem_have_rcpp))
    return(.fsem_kernel_env$fsem_block_xtSx(F, Si, as.integer(nb)))
  n <- nrow(F) / nb
  p <- ncol(F)
  out <- matrix(0, p, p)
  for (k in seq_len(n)) {
    idx <- ((k - 1L) * nb + 1L):((k - 1L) * nb + nb)
    Fk  <- F[idx, , drop = FALSE]
    out <- out + crossprod(Fk, Si %*% Fk)
  }
  out
}

## t(F) %*% kron(I_n, Si) %*% z   (z: length n*nb)
.fsem_xtSz <- function(F, Si, z, nb) {
  if (isTRUE(.fsem_have_rcpp))
    return(.fsem_kernel_env$fsem_block_xtSz(F, Si, as.numeric(z), as.integer(nb)))
  n <- nrow(F) / nb
  p <- ncol(F)
  out <- numeric(p)
  for (k in seq_len(n)) {
    idx <- ((k - 1L) * nb + 1L):((k - 1L) * nb + nb)
    Fk  <- F[idx, , drop = FALSE]
    out <- out + crossprod(Fk, Si %*% z[idx])
  }
  out
}

## ---- closed-form constant penalty blocks ---------------------------------
## Historical penalty P + Ps + Pt (constant given basis + n.b).
.fsem_hist_penalty <- function(basis, n.b) {
  times <- seq(0, 1, length.out = 100)
  E  <- t(eval.basis(times, basis))
  E1 <- t(eval.basis(times, basis, Lfdobj = 1))
  E2 <- t(eval.basis(times, basis, Lfdobj = 2))
  Tn <- length(times)
  G0 <- tcrossprod(E)
  G1 <- tcrossprod(E1)
  G2 <- tcrossprod(E2)
  (kronecker(G1, G1) + kronecker(G2, G0) + kronecker(G0, G2)) / (Tn^2)
}

## Smooth penalty kron(E2 E2^T, sm sm^T) / T  (constant given basis, n.b, sm).
.fsem_smooth_penalty <- function(basis, n.b, sm) {
  times <- seq(0, 1, length.out = 100)
  E2 <- t(eval.basis(times, basis, Lfdobj = 2))
  Tn <- length(times)
  kronecker(tcrossprod(E2), tcrossprod(sm)) / Tn
}

## Precompute pen.fac (per factor model) and pen.sem (per regression) once.
.fsem_precompute_penalties <- function(model, basis, n.b, x.data,
                                       range.min, range.max) {
  bp        <- getbasispenalty(basis)
  hist.pen  <- .fsem_hist_penalty(basis, n.b)

  ## ---- factor model penalties ----
  no.fac <- length(model$mod$factorModel)
  pen.fac.list <- vector("list", no.fac)
  for (j in seq_len(no.fac)) {
    pen.fac <- bp
    for (jj in seq_along(model$mod$factorModel[[j]]$effect)) {
      eff <- model$mod$factorModel[[j]]$effect[jj]
      pe <- NULL
      if (eff %in% c("fixed_concurrent", "fixed_historical")) pe <- NULL
      if (eff == "concurrent")  pe <- bp
      if (eff == "historical")  pe <- hist.pen
      if (is.null(pe)) next
      pen.fac <- Matrix::bdiag(pen.fac, pe)
    }
    pen.fac.list[[j]] <- as.matrix(pen.fac)
  }

  ## ---- regression (SEM) penalties ----
  no.reg <- length(model$mod$regression)
  pen.sem.list <- vector("list", no.reg)
  for (j in seq_len(no.reg)) {
    effs <- model$mod$regression[[j]]$effect
    if (length(effs) == 0) { pen.sem.list[[j]] <- NULL; next }
    pen.sem <- 0
    for (jj in seq_along(effs)) {
      eff <- effs[jj]
      pe <- NULL
      if (eff == "concurrent") pe <- bp
      if (eff == "historical") pe <- hist.pen
      if (eff == "smooth") {
        cov.name <- model$mod$regression[[j]]$covariate[jj]
        basis3 <- create.bspline.basis(rangeval = range(x.data[[cov.name]]), nbasis = n.b)
        sm <- t(eval.basis(x.data[[cov.name]], basis3))
        pe <- .fsem_smooth_penalty(basis, n.b, sm)
      }
      if (eff == "linear") pe <- bp
      if (is.null(pe)) next
      pen.sem <- Matrix::bdiag(pen.sem, pe)
    }
    pen.sem <- pen.sem[-1, -1]
    pen.sem.list[[j]] <- as.matrix(pen.sem)
  }

  list(fac = pen.fac.list, sem = pen.sem.list)
}

######################################################################
## em.estimation.optimised : same signature & algorithm as em.estimation
######################################################################
em.estimation.optimised <- function(model, data, x.data = NULL,
                                     initial.parameter, n.b, n.em = 100, n.monte = 100,
                                     s.p = c("min", "sd"), range.min = NULL, range.max = NULL,
                                     design = c("regular", "irregular", "regular.truncated", "regular.missing"),
                                     plot.progress = FALSE,
                                     save.every = NULL) {

  s.p = match.arg(s.p)
  if (is.null(range.min)) range.min = 0
  if (is.null(range.max)) range.max = 1
  if (is.null(s.p)) s.p = "min"

  .fsem_load_kernels()

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

  ## ---- OPTIMISATION 1: precompute all constant penalty matrices ONCE ----
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

  ###############Algorithm#####################################
  for (i in 1:no.em) {

    eta.distt <- eta.distribution(model, param = parameterr, n.b, range.min, range.max, data, x.data, omegaMatrix = omegaMatrix)
    paramz1 <- params.z1(model, param = parameterr, n.b, range.min, range.max, eta.dist = eta.distt, data, omegaMatrix = omegaMatrix, design = design)
    time.start <- Sys.time()

    ################fac parameter estimations#######################
    generators2.z <- list()
    generators2.w <- list()
    for (i1 in 1:samples) {
      samp <- i1
      gg <- list()
      mu.gg <- list()
      generators2.z[[paste0(samp)]] <- list()
      generators2.w[[paste0(samp)]] <- list()
      zz <- 0
      w.tot <- data$.value[data$.id == samp]
      paramz <- params.z(model, param = paramz1, n.b = n.b,
                         range.min = range.min, range.max = range.max,
                         sample.no = i1, w = w.tot, data, eval.t = eval.t)
      paramz11 <- list()
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
        ## NOTE: kept arithmetically identical to the original so that the
        ## eigen-decomposition inside mvrnorm() (sign of eigenvectors) is
        ## bit-identical and RNG realisations match exactly.
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
        mean <- rowMeans(sapply(1:samples, function(k) { eta.l[[m]][[k]][[i2]] }))
        for (k in 1:samples) eta.l[[m]][[k]][[i2]] <- eta.l[[m]][[k]][[i2]] - mean
      }
    }

    if (design.regular == TRUE) {
      soll <- ginv(eval.t$eval[[1]][[1]] %*% t(eval.t$eval[[1]][[1]])) %*% eval.t$eval[[1]][[1]]
    }
    delta.fac <- c()
    hist.fac <- list()
    for (j in 1:no.fac) {
      aa <- c()
      for (k in 1:dim(diagMatrix2[[paste0("f", j)]])[1]) {
        aa[k] <- if (all(diagMatrix2[[paste0("f", j)]][k, ] == 0)) 0 else 1
      }
      ## OPTIMISATION 1: use precomputed penalty
      pen.fac <- pen.fac.list[[j]]

      f <- list(); f2 <- list()
      for (m in 1:mont) {
        f[[m]] <- list(); f2[[m]] <- list()
        for (k in 1:samples) {
          if (design.regular == FALSE) {
            soll <- ginv(eval.t$eval[[k]][[j]] %*% t(eval.t$eval[[k]][[j]])) %*% eval.t$eval[[k]][[j]]
          }
          f[[m]][[paste0(k)]] <- diag(n.b)
          for (i2 in 1:no.lat) {
            et <- kronecker(diag(eval.t$time.no[[k]][[j]]), t(eta.l[[m]][[k]][[i2]]))
            fac <- which(model$mod$factorModel[[j]]$factor == model$var$latents[i2])
            if (!length(fac) == 0) {
              if (model$mod$factorModel[[j]]$effect[fac] == "concurrent") {
                ef <- soll %*% (et %*% omegaMatrix$omega1[[k]][[j]]); f2[[m]][[paste0(k)]] <- 0
              }
              if (model$mod$factorModel[[j]]$effect[fac] == "historical") {
                ef <- soll %*% (et %*% omegaMatrix$omega2[[k]][[j]]); f2[[m]][[paste0(k)]] <- 0
              }
              if (model$mod$factorModel[[j]]$effect[fac] == "fixed_concurrent") {
                ef <- NULL; f2[[m]][[paste0(k)]] <- eta.l[[m]][[k]][[i2]]
              }
              if (model$mod$factorModel[[j]]$effect[fac] == "fixed_historical") {
                ef <- NULL; f2[[m]][[paste0(k)]] <- soll %*% omegaMatrix$delta[[k]][[j]] %*% eta.l[[m]][[k]][[i2]]
              }
            } else {
              ef <- NULL
            }
            f[[m]][[paste0(k)]] <- cbind(f[[m]][[paste0(k)]], ef)
          }
        }
      }

      if (i > 0) {
        delt_sta <- c(0.001, 0.01, 0.05, 0.1) * log(samples)
        delt <- delt_sta
        samp.train <- sort(sample(1:samples, round(0.7 * samples)))
        samp.test  <- sort((1:samples)[-samp.train])
        sig <- parameterr$sigma.fac[[j]]
        ## OPTIMISATION: build stacked matrices without rbind-growth
        ff <- list(); zz <- list(); ff2 <- list()
        for (m in 1:mont) {
          ff[[m]]  <- do.call(rbind, lapply(samp.train, function(i1) f[[m]][[paste0(i1)]]))
          zz[[m]]  <- unlist(lapply(samp.train, function(i1) generators2.z[[paste0(i1)]][[j]][[m]]))
          ff2[[m]] <- unlist(lapply(samp.train, function(i1) f2[[m]][[paste0(i1)]]))
        }
        ff.test <- list(); zz.test <- list(); ff2.test <- list()
        for (m in 1:mont) {
          ff.test[[m]]  <- do.call(rbind, lapply(samp.test, function(i1) f[[m]][[paste0(i1)]]))
          zz.test[[m]]  <- unlist(lapply(samp.test, function(i1) generators2.z[[paste0(i1)]][[j]][[m]]))
          ff2.test[[m]] <- unlist(lapply(samp.test, function(i1) f2[[m]][[paste0(i1)]]))
        }
        per <- c()
        A.test <- vector("list", mont); b.test <- vector("list", mont)
        for (l in 1:length(delt)) {
          for (it in 1:no.it) {
            ## OPTIMISATION 2: block-diagonal inverse (invert only nb x nb sig)
            Si <- solve(sig)
            lam.sum <- 0; lam1.sum <- 0
            for (m in 1:mont) {
              lam  <- .fsem_xtSx(ff[[m]], Si, n.b) + delt[l] * pen.fac
              lam1 <- .fsem_xtSz(ff[[m]], Si, as.vector(zz[[m]]) - ff2[[m]], n.b)
              lam.sum  <- lam.sum + 1 / mont * lam
              lam1.sum <- lam1.sum + 1 / mont * lam1
            }
            lambda1 <- solve(lam.sum) %*% lam1.sum
            lambda2.1 <- lambda1
            sigma1 <- 0
            for (i1 in samp.train) {
              ss2 <- 0
              for (m in 1:mont) {
                ss <- as.vector(generators2.z[[paste0(i1)]][[j]][[m]]) - f2[[m]][[paste0(i1)]] - f[[m]][[paste0(i1)]] %*% lambda2.1
                ss2 <- ss2 + 1 / mont * (ss %*% t(ss))
              }
              sigma1 <- sigma1 + 1 / length(samp.train) * ss2
            }
            sig <- sigma1
          }
          Si.test <- solve(sig)
          likeli <- 0
          for (m in 1:mont) {
            A.test[[m]] <- .fsem_xtSx(ff.test[[m]], Si.test, n.b)
            b.test[[m]] <- .fsem_xtSz(ff.test[[m]], Si.test, as.vector(zz.test[[m]]) - ff2.test[[m]], n.b)
            vall <- (-t(lambda2.1) %*% b.test[[m]]) +
              t(lambda2.1) %*% (A.test[[m]] + delt[l] * pen.fac) %*% lambda2.1
            likeli <- likeli + vall / mont
          }
          per[l] <- as.vector(likeli)
        }
        delt.min <- delt[which(per == min(per))]
        ## A.test/b.test currently hold values from the last delt (matches original 'sol')
        ex1 <- c()
        for (m in 1:mont) {
          vall <- (-t(lambda2.1) %*% b.test[[m]]) +
            t(lambda2.1) %*% (A.test[[m]] + delt.min * pen.fac) %*% lambda2.1
          ex1[m] <- as.vector(vall)
        }
        bound <- mean(ex1) + sd(ex1)
        per <- c()
        for (l in 1:length(delt)) {
          likeli <- 0
          for (m in 1:mont) {
            vall <- (-t(lambda2.1) %*% b.test[[m]]) +
              t(lambda2.1) %*% (A.test[[m]] + delt[l] * pen.fac) %*% lambda2.1
            likeli <- likeli + 1 / mont * as.vector(vall)
          }
          per[l] <- (likeli <= bound)
        }
        d.t <- delt[which(per)]
        deltaa <- if (s.p == "min") delt.min else max(d.t)
      } else {
        deltaa <- ifelse(n.b > 2, 0.01, 0)
      }
      delta.fac[j] <- deltaa
      ff <- list(); zz <- list(); ff2 <- list()
      for (m in 1:mont) {
        ff[[m]]  <- do.call(rbind, lapply(1:samples, function(i1) f[[m]][[paste0(i1)]]))
        zz[[m]]  <- unlist(lapply(1:samples, function(i1) generators2.z[[paste0(i1)]][[j]][[m]]))
        ff2[[m]] <- unlist(lapply(1:samples, function(i1) f2[[m]][[paste0(i1)]]))
      }
      sig <- parameterr$sigma.fac[[j]]
      for (it in 1:no.it) {
        Si <- solve(sig)
        lam.sum <- 0; lam1.sum <- 0
        for (m in 1:mont) {
          lam  <- .fsem_xtSx(ff[[m]], Si, n.b) + deltaa * pen.fac
          lam1 <- .fsem_xtSz(ff[[m]], Si, as.vector(zz[[m]]) - ff2[[m]], n.b)
          lam.sum  <- lam.sum + 1 / mont * lam
          lam1.sum <- lam1.sum + 1 / mont * lam1
        }
        lambda1 <- solve(lam.sum) %*% lam1.sum
        lambda2.1 <- lambda1
        sigma1 <- 0
        for (i1 in 1:samples) {
          ss2 <- 0
          for (m in 1:mont) {
            ss <- as.vector(generators2.z[[paste0(i1)]][[j]][[m]]) - f2[[m]][[paste0(i1)]] - f[[m]][[paste0(i1)]] %*% lambda2.1
            ss2 <- ss2 + 1 / mont * (ss %*% t(ss))
          }
          sigma1 <- sigma1 + 1 / samples * ss2
        }
        sig <- sigma1
      }

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

    ################sem parameter estimations##################
    no.reg <- length(model$mod$regression)
    lat2 <- model$var$latents
    obs <- model$var$observed
    no.lat2 <- length(lat2)
    no.obs2 <- length(obs)
    soll <- ginv(evall2 %*% t(evall2)) %*% evall2
    hist.sem <- list()
    for (j in 1:no.reg) {
      ## Reset cross-iteration flags: et.ex/x.ex are set (via <-) only when a
      ## regression has no latent / no observed covariate. Because exists()
      ## checks the whole function scope, a stale flag from a previous
      ## regression would wrongly strip a design column in cbind() and make the
      ## design Gram matrix non-conformable with pen.sem. Clear them each pass.
      if (exists("et.ex")) rm(et.ex)
      if (exists("x.ex"))  rm(x.ex)
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
      ## OPTIMISATION 1: use precomputed SEM penalty (built only if needed)
      if (length(which(aa == 1)) > 0) {
        pen.sem <- pen.sem.list[[j]]
      }
      g.eta <- NULL; g.x <- NULL
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
          et.ex <- "NULL"
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
            g.x <- lapply(1:samples, function(i1) (g.x[[i1]] = 0)); x.ex <- "NULL"
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
          g <- list(); gg <- list(); ett <- list()
          if (!is.null(g.eta) | !is.null(g.x)) {
            for (m in 1:mont) {
              g[[m]] <- list(); gg[[m]] <- 0; ett[[m]] <- 0
              for (i1 in samp.train) {
                g[[m]][[i1]] <- cbind(g.eta[[m]][[i1]], g.x[[i1]])
                if (exists("et.ex")) g[[m]][[i1]] <- g[[m]][[i1]][, -1]
                if (exists("x.ex")) g[[m]][[i1]] <- g[[m]][[i1]][, -ncol(g[[m]][[i1]])]
                gg[[m]] <- rbind(gg[[m]], g[[m]][[i1]])
                ett[[m]] <- c(ett[[m]], eta.l[[m]][[i1]][[lat.count]])
              }
              gg[[m]] <- gg[[m]][-1, ]; ett[[m]] <- ett[[m]][-1]
            }
            g.test <- list(); gg.test <- list(); ett.test <- list()
            for (m in 1:mont) {
              g.test[[m]] <- list(); gg.test[[m]] <- 0; ett.test[[m]] <- 0
              for (i1 in samp.test) {
                g.test[[m]][[i1]] <- cbind(g.eta[[m]][[i1]], g.x[[i1]])
                if (exists("et.ex")) g.test[[m]][[i1]] <- g.test[[m]][[i1]][, -1]
                if (exists("x.ex")) g.test[[m]][[i1]] <- g.test[[m]][[i1]][, -ncol(g.test[[m]][[i1]])]
                gg.test[[m]] <- rbind(gg.test[[m]], g.test[[m]][[i1]])
                ett.test[[m]] <- c(ett.test[[m]], eta.l[[m]][[i1]][[lat.count]])
              }
              gg.test[[m]] <- gg.test[[m]][-1, ]; ett.test[[m]] <- ett.test[[m]][-1]
            }
          }

          per <- c()
          A.test <- vector("list", mont); b.test <- vector("list", mont)
          for (l in 1:length(delt)) {
            for (it in 1:no.it) {
              Si <- solve(sig)
              gam.sum <- 0; gam1.sum <- 0
              for (m in 1:mont) {
                gam  <- .fsem_xtSx(gg[[m]], Si, n.b) + delt[l] * pen.sem
                gam1 <- .fsem_xtSz(gg[[m]], Si, ett[[m]], n.b)
                gam.sum  <- gam.sum + 1 / mont * gam
                gam1.sum <- gam1.sum + 1 / mont * gam1
              }
              gamma1 <- solve(gam.sum) %*% gam1.sum
              gamma2.1 <- gamma1
              sigma1 <- 0
              for (i1 in samp.train) {
                ss2 <- 0
                for (m in 1:mont) {
                  if (is.null(g.eta) & is.null(g.x)) {
                    ss <- eta.l[[m]][[i1]][[lat.count]]
                  } else {
                    ss <- eta.l[[m]][[i1]][[lat.count]] - g[[m]][[i1]] %*% gamma2.1
                  }
                  ss2 <- ss2 + 1 / mont * ss %*% t(ss)
                }
                sigma1 <- sigma1 + 1 / length(samp.train) * ss2
              }
              sig <- sigma1
            }
            Si.test <- solve(sig)
            val <- 0
            for (m in 1:mont) {
              A.test[[m]] <- .fsem_xtSx(gg.test[[m]], Si.test, n.b)
              b.test[[m]] <- .fsem_xtSz(gg.test[[m]], Si.test, as.vector(ett.test[[m]]), n.b)
              vall <- (-t(gamma2.1) %*% b.test[[m]]) +
                t(gamma2.1) %*% (A.test[[m]] + delt[l] * pen.sem) %*% gamma2.1
              val <- val + 1 / mont * vall
            }
            per[l] <- as.vector(val)
          }
          delt.min <- delt[which(per == min(per))]
          ex1 <- c()
          for (m in 1:mont) {
            vall <- (-t(gamma2.1) %*% b.test[[m]]) +
              t(gamma2.1) %*% (A.test[[m]] + delt.min * pen.sem) %*% gamma2.1
            ex1[m] <- as.vector(vall)
          }
          bound <- mean(ex1) + sd(ex1)
          per <- c()
          for (l in 1:length(delt)) {
            likeli <- 0
            for (m in 1:mont) {
              vall <- (-t(gamma2.1) %*% b.test[[m]]) +
                t(gamma2.1) %*% (A.test[[m]] + delt[l] * pen.sem) %*% gamma2.1
              likeli <- likeli + 1 / mont * as.vector(vall)
            }
            per[l] <- (likeli <= bound)
          }
          d.t <- delt[which(per)]
          deltaa <- if (s.p == "min") delt.min else max(d.t)
        } else {
          deltaa <- 0
        }
        delta.sem <- deltaa
        g <- list(); gg <- list(); ett <- list()
        if (!is.null(g.eta) | !is.null(g.x))
          for (m in 1:mont) {
            g[[m]] <- list(); gg[[m]] <- 0; ett[[m]] <- 0
            for (i1 in 1:samples) {
              g[[m]][[i1]] <- cbind(g.eta[[m]][[i1]], g.x[[i1]])
              if (exists("et.ex")) g[[m]][[i1]] <- g[[m]][[i1]][, -1]
              if (exists("x.ex")) g[[m]][[i1]] <- g[[m]][[i1]][, -ncol(g[[m]][[i1]])]
              gg[[m]] <- rbind(gg[[m]], g[[m]][[i1]])
              ett[[m]] <- c(ett[[m]], eta.l[[m]][[i1]][[lat.count]])
            }
            gg[[m]] <- gg[[m]][-1, ]; ett[[m]] <- ett[[m]][-1]
          }
        sig <- parameterr$sigma.sem[[lat.count]]
        for (it in 1:no.it) {
          Si <- solve(sig)
          gam.sum <- 0; gam1.sum <- 0
          for (m in 1:mont) {
            gam  <- .fsem_xtSx(gg[[m]], Si, n.b) + deltaa * pen.sem
            gam1 <- .fsem_xtSz(gg[[m]], Si, ett[[m]], n.b)
            gam.sum  <- gam.sum + 1 / mont * gam
            gam1.sum <- gam1.sum + 1 / mont * gam1
          }
          gamma1 <- solve(gam.sum) %*% gam1.sum
          gamma2.1 <- gamma1
          sigma1 <- 0
          for (i1 in 1:samples) {
            ss2 <- 0
            for (m in 1:mont) {
              if (is.null(g.eta) & is.null(g.x)) {
                ss <- eta.l[[m]][[i1]][[lat.count]]
              } else {
                ss <- eta.l[[m]][[i1]][[lat.count]] - g[[m]][[i1]] %*% gamma2.1
              }
              ss2 <- ss2 + 1 / mont * ss %*% t(ss)
            }
            sigma1 <- sigma1 + 1 / samples * ss2
          }
          sig <- sigma1
        }
      } else {
        sigma1 <- 0
        for (i1 in 1:samples) {
          ss2 <- 0
          for (m in 1:mont) {
            if (is.null(g.eta) & is.null(g.x)) {
              ss <- eta.l[[m]][[i1]][[lat.count]]
            } else {
              ss <- eta.l[[m]][[i1]][[lat.count]] - g[[m]][[i1]] %*% gamma2.1
            }
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
        hist.sem[[lat.count]] <- list(g.eta = g.eta, g.x = g.x, g = g, deltaa = deltaa, pen.sem = pen.sem)
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

    ###convergence checking##################
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
      if (plot.progress) {
        plot(1:length(toler_mean_sig[-1]), toler_mean_sig[-1], "l")
        plot(1:length(toler_mean[-1]), toler_mean[-1], "l")
      }
    }
    #####################################

    if (plot.progress) {
      t <- seq(0, 1, length.out = 200)
      par(mfrow = c(2, 2))
      plot(t, params.eval$coef.fac.std$z1$eta, col = "red", "l")
      plot(t, params.eval$coef.fac.std$z2$eta, col = "red", "l")
      plot(t, params.eval$coef.fac.std$z3$eta, col = "red", "l")
      plot(t, params.eval$coef.fac.std$z4$eta, col = "red", "l")
      plot(t, params.eval$coef.fac.std$z5$eta, col = "red", "l")
      par(mfrow = c(2, 2))
      plot(t, params.eval$coef.sem.std$eta$gender, col = "red", "l")
      plot(t, params.eval$coef.sem.std$eta$cancer, col = "red", "l")
      plot(t, params.eval$coef.sem.std$eta$genCancer, col = "red", "l")
      par(mfrow = c(2, 2))
      plot(t, params.eval$sigma.sem.eigvec.std$eta[, 1], col = "red", "l")
      plot(t, params.eval$sigma.sem.eigvec.std$eta[, 2], col = "red", "l")
      plot(t, params.eval$sigma.sem.eigvec.std$eta[, 3], col = "red", "l")
    }

    ##########################################################################
    param_prev <- params.eval
    sig_prev <- parameterr$sigma.error
    res.k = res.k + 1
    res.list[[res.k]] = list(parameters = parameters, parameters.eval = params.eval, timing = c(time.start = time.start, time.end = time.end))

    if (plot.progress) print(parameterr$sigma.error)

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
