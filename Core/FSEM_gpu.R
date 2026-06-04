######################################################################
## FSEM_gpu.R
##
## GPU re-implementation of em.estimation.optimised.revised() using the
## `torch` package (libtorch + CUDA).  Defines:
##
##     em.estimation.optimised.gpu(...)
##
## with the SAME signature as em.estimation.optimised.revised() (it also
## applies the gamma jitter that fixes the coupled-latent EM degeneracy),
## and the SAME statistical algorithm as em.estimation.optimised().
##
## WHAT RUNS ON THE GPU
## --------------------
## The compute-dominant part of every EM iteration is the M-step:
## penalised weighted least squares (IRLS) with a 4-point delta cross
## validation, repeated for every factor model and every structural
## regression, with an inner reduction over the Monte-Carlo draws and the
## sample blocks.  This is the documented hot path (see FSEM_optimised.R).
##
## In this implementation the ENTIRE M-step numeric engine is performed on
## the GPU with `torch` tensors that stay RESIDENT on the device across all
## IRLS / CV iterations for a given block:
##   * the stacked per-(monte, sample) design tensors are uploaded ONCE per
##     factor / regression,
##   * the block quadratic / bilinear forms  t(F) kron(I,Si) F  and
##     t(F) kron(I,Si) z  are evaluated as single batched `torch_einsum`
##     reductions over (monte x sample x basis),
##   * the normal-equation solves, the sigma updates, the delta-CV
##     likelihood evaluations and the bound/selection arithmetic all run on
##     the device.  No host<->device transfer happens inside the IRLS loop.
##
## WHAT STAYS ON THE CPU
## ---------------------
## The model-structure bookkeeping that is NOT a numeric hot spot is reused
## verbatim from the protected helpers in Core/FSEM.R (eta.distribution,
## params.z1, params.z, generator.z, ex.estimations,
## parameter.estimated.evaluted) and Core/FSEM_optimised.R (penalty
## precomputation).  These build the per-sample posterior moments, draw the
## Monte-Carlo latent samples and assemble the (model-coupled) design
## matrices.  A literal "100% on GPU" port is not possible without rewriting
## those protected files; everything that is arithmetic-heavy is on the GPU.
##
## FIDELITY
## --------
## Results are statistically equivalent to em.estimation.optimised.revised
## but NOT bit-identical: GPU floating-point reduction order and (CPU) RNG
## differ.  Use the same seed for the gamma jitter to reproduce the GPU run.
##
## DEVICE FALLBACK
## ---------------
## If torch/CUDA is unavailable (or CUDA kernels cannot execute on the
## installed GPU), the code automatically falls back to CPU torch tensors so
## the function still runs and returns valid estimates.
##
## Requirements: source AFTER Core/FSEM.R and Core/FSEM_optimised.R.
######################################################################

## ---- GPU device / dtype state -------------------------------------------
.fsem_gpu_state <- new.env(parent = globalenv())
.fsem_gpu_state$device <- NULL   # torch_device, set by .fsem_gpu_init
.fsem_gpu_state$dtype  <- NULL
.fsem_gpu_state$desc   <- NULL

## Initialise torch and choose the compute device.
##   prefer = "cuda"  : try CUDA, fall back to CPU if a real kernel fails.
##   prefer = "cpu"   : force CPU tensors.
## A small matmul is actually executed on CUDA to confirm that kernels run
## on the installed GPU (cuda_is_available() alone is not sufficient on very
## new architectures), otherwise we fall back to CPU.
.fsem_gpu_init <- function(prefer = c("cuda", "cpu"), verbose = TRUE) {
  prefer <- match.arg(prefer)
  if (!requireNamespace("torch", quietly = TRUE))
    stop("em.estimation.optimised.gpu requires the 'torch' package. Install it with install.packages('torch') and torch::install_torch().")
  if (!torch::torch_is_installed())
    stop("The torch C++ backend (libtorch) is not installed. Run torch::install_torch() (set CUDA=\"12.8\" for an RTX 50-series GPU).")

  dtype <- torch::torch_float64()
  dev   <- torch::torch_device("cpu")
  desc  <- "cpu"

  if (prefer == "cuda" && isTRUE(try(torch::cuda_is_available(), silent = TRUE))) {
    ok <- tryCatch({
      d  <- torch::torch_device("cuda")
      a  <- torch::torch_randn(64, 64, dtype = dtype, device = d)
      b  <- torch::torch_mm(a, a)
      v  <- as.numeric(b[1, 1]$cpu())
      is.finite(v)
    }, error = function(e) {
      if (verbose) message("  [gpu] CUDA present but kernel test failed (",
                           conditionMessage(e), "); using CPU torch tensors.")
      FALSE
    })
    if (isTRUE(ok)) { dev <- torch::torch_device("cuda"); desc <- "cuda" }
  }

  .fsem_gpu_state$device <- dev
  .fsem_gpu_state$dtype  <- dtype
  .fsem_gpu_state$desc   <- desc
  if (verbose) message("  [gpu] em.estimation.optimised.gpu compute device: ", desc)
  invisible(desc)
}

## R matrix / vector  ->  torch tensor on the active device (float64)
.fsem_t <- function(x) {
  torch::torch_tensor(x, dtype = .fsem_gpu_state$dtype, device = .fsem_gpu_state$device)
}

## ---------------------------------------------------------------------
## .fsem_gpu_penalized_cv_fit
## ---------------------------------------------------------------------
## GPU engine for the penalised-WLS + 4-point delta cross-validation that
## the original performs for every factor model and every structural
## regression.  It reproduces, on the device, exactly the arithmetic of the
## CPU blocks in em.estimation.optimised() (factor: lines building
## ff/zz/ff2 + IRLS + CV + final; SEM: the analogous gg/ett block).
##
## Inputs (all built CPU-side, model-coupled, then uploaded ONCE):
##   D : list over m (1..mont) of list over k (1..samples) of (nb x p) design
##   R : list over m of list over k of length-nb response vector
##   O : list over m of list over k of length-nb offset (or scalar 0)
##   pen        : (p x p) penalty matrix
##   sig.init   : (nb x nb) initial residual covariance
##   mont, samples, nb : dimensions
##   no.it      : IRLS iterations per (delta) / for the final fit
##   s.p        : "min" or "sd"  delta-selection rule
##   delt       : numeric vector of candidate penalties (length L)
##   samp.train, samp.test : integer index vectors into 1..samples
##   do.cv      : if FALSE, skip CV and fit with delt[?]=0 (matches i==0 path)
##
## Returns list(lambda = numeric(p), sig = matrix(nb,nb), deltaa = numeric).
## All heavy work is on the GPU; only tiny scalars are copied back to host.
.fsem_gpu_penalized_cv_fit <- function(D, R, O, pen, sig.init,
                                       mont, samples, nb, no.it,
                                       s.p, delt, samp.train, samp.test,
                                       do.cv = TRUE) {
  p   <- ncol(D[[1]][[1]])
  dev <- .fsem_gpu_state$device
  dty <- .fsem_gpu_state$dtype

  pen.t  <- .fsem_t(pen)
  sig0.t <- .fsem_t(sig.init)

  ## --- assemble a (mont, S, nb, p) design tensor + (mont, S, nb) target ---
  ## target t = R - O  (offset of scalar 0 contributes nothing, exactly as
  ## the original recycling of length-`samples` zero vectors).
  build_tensors <- function(idx) {
    S <- length(idx)
    Darr <- array(0, dim = c(mont, S, nb, p))
    Tarr <- array(0, dim = c(mont, S, nb))
    for (m in seq_len(mont)) {
      Dm <- D[[m]]; Rm <- R[[m]]; Om <- O[[m]]
      for (s in seq_len(S)) {
        k <- idx[s]
        Darr[m, s, , ] <- Dm[[k]]
        rk <- as.numeric(Rm[[k]])
        ok <- Om[[k]]
        if (length(ok) == nb) rk <- rk - as.numeric(ok)  ## scalar 0 -> no-op
        Tarr[m, s, ] <- rk
      }
    }
    list(D = .fsem_t(Darr), T = .fsem_t(Tarr), S = S)
  }

  ## --- one penalised-WLS IRLS fit (no.it iterations) on a resident set ---
  ## Dt: (mont,S,nb,p)  Tt: (mont,S,nb)  delta: scalar  sig.t: (nb,nb) start
  ## sig.t is UPDATED in place across iterations (returned).  Returns lambda.
  irls <- function(Dt, Tt, S, delta, sig.t) {
    lambda <- NULL
    dpen <- torch::torch_mul(pen.t, delta)
    for (it in seq_len(no.it)) {
      Si <- torch::linalg_inv(sig.t)                                   # (nb,nb)
      ## A = (1/mont) sum_{m,k} D^T Si D  + delta*pen
      A <- torch::torch_einsum("mkap,ab,mkbq->pq", list(Dt, Si, Dt))
      A <- torch::torch_add(torch::torch_div(A, mont), dpen)
      ## b = (1/mont) sum_{m,k} D^T Si t
      b <- torch::torch_einsum("mkap,ab,mkb->p", list(Dt, Si, Tt))
      b <- torch::torch_div(b, mont)
      lambda <- torch::linalg_solve(A, b$unsqueeze(2))$squeeze(2)      # (p)
      ## residual + sigma update: sig = (1/(mont*S)) sum r r^T
      pred <- torch::torch_einsum("mkap,p->mka", list(Dt, lambda))
      r    <- torch::torch_sub(Tt, pred)
      sig.t <- torch::torch_div(
                 torch::torch_einsum("mka,mkb->ab", list(r, r)),
                 mont * S)
    }
    list(lambda = lambda, sig = sig.t)
  }

  ## test-set per-monte A_m (mont,p,p) and b_m (mont,p) for a given Si
  testAb <- function(Dt, Tt, Si) {
    A <- torch::torch_einsum("mkap,ab,mkbq->mpq", list(Dt, Si, Dt))
    b <- torch::torch_einsum("mkap,ab,mkb->mp",  list(Dt, Si, Tt))
    list(A = A, b = b)
  }
  ## mean_m ( -b_m . lambda + lambda^T (A_m + d*pen) lambda )
  cv_obj <- function(A, b, lambda, d) {
    t1 <- torch::torch_einsum("mp,p->m", list(b, lambda))
    Ap <- torch::torch_add(A, torch::torch_mul(pen.t, d))
    t2 <- torch::torch_einsum("p,mpq,q->m", list(lambda, Ap, lambda))
    torch::torch_sub(t2, t1)                                           # (mont)
  }

  full <- build_tensors(seq_len(samples))

  if (!do.cv) {
    fit <- irls(full$D, full$T, full$S, 0, sig0.t$clone())
    return(list(lambda = as.numeric(fit$lambda$cpu()),
                sig    = as.matrix(fit$sig$cpu()),
                deltaa = 0))
  }

  tr <- build_tensors(samp.train)
  te <- build_tensors(samp.test)
  L  <- length(delt)

  per      <- numeric(L)
  sig.t    <- sig0.t$clone()          # sig carries across delta candidates
  lambda.l <- NULL; At <- NULL; bt <- NULL
  for (l in seq_len(L)) {
    fit    <- irls(tr$D, tr$T, tr$S, delt[l], sig.t)
    sig.t  <- fit$sig
    lambda.l <- fit$lambda
    Si.te  <- torch::linalg_inv(sig.t)
    ab     <- testAb(te$D, te$T, Si.te)
    At <- ab$A; bt <- ab$b
    per[l] <- as.numeric(torch::torch_mean(cv_obj(At, bt, lambda.l, delt[l]))$cpu())
  }
  delt.min <- delt[which(per == min(per))][1]

  ## bound = mean+sd of the per-monte objective at delt.min (using last At/bt)
  ex1   <- as.numeric(cv_obj(At, bt, lambda.l, delt.min)$cpu())
  bound <- mean(ex1) + sd(ex1)
  per2  <- logical(L)
  for (l in seq_len(L)) {
    val <- as.numeric(torch::torch_mean(cv_obj(At, bt, lambda.l, delt[l]))$cpu())
    per2[l] <- (val <= bound)
  }
  d.t    <- delt[per2]
  deltaa <- if (s.p == "min") delt.min else max(d.t)

  ## final fit on ALL samples, sigma reset to the same initial value
  fit <- irls(full$D, full$T, full$S, deltaa, sig0.t$clone())
  list(lambda = as.numeric(fit$lambda$cpu()),
       sig    = as.matrix(fit$sig$cpu()),
       deltaa = deltaa)
}

######################################################################
## em.estimation.optimised.gpu
##   - same signature as em.estimation.optimised.revised
##   - same gamma jitter (fixes coupled-latent degeneracy)
##   - E-step + design assembly via the protected CPU helpers
##   - M-step penalised-WLS / CV / sigma updates on the GPU
######################################################################
em.estimation.optimised.gpu <- function(model, data, x.data = NULL,
                                         initial.parameter, n.b,
                                         n.em = 100, n.monte = 100,
                                         s.p = c("min", "sd"),
                                         range.min = NULL, range.max = NULL,
                                         design = c("regular", "irregular",
                                                    "regular.truncated", "regular.missing"),
                                         plot.progress = FALSE,
                                         save.every = NULL,
                                         gamma.init.sd = 0.3,
                                         seed = NULL,
                                         gpu.device = c("cuda", "cpu"),
                                         gpu.verbose = TRUE) {

  s.p = match.arg(s.p)
  design = match.arg(design)
  gpu.device = match.arg(gpu.device)
  if (is.null(range.min)) range.min = 0
  if (is.null(range.max)) range.max = 1

  ## ---- initialise the GPU compute device ----
  .fsem_gpu_init(prefer = gpu.device, verbose = gpu.verbose)

  ## ---- gamma jitter (identical to em.estimation.optimised.revised) ----
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
      if (length(g) > 0)
        init.perturbed$gamma.param[[k]] <- g + rnorm(length(g), 0, gamma.init.sd)
    }
  }

  .fsem_load_kernels()  ## only used by the (unused here) CPU block fallback

  omegaMatrix <- omega.matrix(n.b, range.min, range.max, data)
  result <- list()
  no.em <- n.em
  no.it <- 10
  parameterr <- init.perturbed
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
  design.regular <- (design == "regular")

  penalties <- .fsem_precompute_penalties(model, basis, n.b, x.data, range.min, range.max)
  pen.fac.list <- penalties$fac
  pen.sem.list <- penalties$sem

  pb <- txtProgressBar(0, no.em, style = 3)
  pbk <- 0

  res.list <- list()
  res.k <- 0
  toler_mean <- c(); toler_mean_sig <- c()

  ###############Algorithm#####################################
  for (i in 1:no.em) {

    eta.distt <- eta.distribution(model, param = parameterr, n.b, range.min, range.max, data, x.data, omegaMatrix = omegaMatrix)
    paramz1 <- params.z1(model, param = parameterr, n.b, range.min, range.max, eta.dist = eta.distt, data, omegaMatrix = omegaMatrix, design = design)
    time.start <- Sys.time()

    ################ E-step: Monte-Carlo latent samples (CPU helpers) ########
    generators2.z <- list()
    generators2.w <- list()
    for (i1 in 1:samples) {
      samp <- i1
      gg <- list(); mu.gg <- list()
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

    ################ M-step: factor parameter estimation (GPU) ##############
    delta.fac <- c()
    hist.fac <- list()
    for (j in 1:no.fac) {
      aa <- c()
      for (k in 1:dim(diagMatrix2[[paste0("f", j)]])[1]) {
        aa[k] <- if (all(diagMatrix2[[paste0("f", j)]][k, ] == 0)) 0 else 1
      }
      pen.fac <- pen.fac.list[[j]]

      ## assemble the (model-coupled) design / response / offset lists
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

      ## response / offset lists indexed [[m]][[k]] (k = 1..samples)
      Dlist <- vector("list", mont); Rlist <- vector("list", mont); Olist <- vector("list", mont)
      for (m in 1:mont) {
        Dlist[[m]] <- vector("list", samples)
        Rlist[[m]] <- vector("list", samples)
        Olist[[m]] <- vector("list", samples)
        for (k in 1:samples) {
          Dlist[[m]][[k]] <- f[[m]][[paste0(k)]]
          Rlist[[m]][[k]] <- as.vector(generators2.z[[paste0(k)]][[j]][[m]])
          Olist[[m]][[k]] <- f2[[m]][[paste0(k)]]
        }
      }

      delt_sta <- c(0.001, 0.01, 0.05, 0.1) * log(samples)
      samp.train <- sort(sample(1:samples, round(0.7 * samples)))
      samp.test  <- sort((1:samples)[-samp.train])

      gpu.fit <- .fsem_gpu_penalized_cv_fit(
        D = Dlist, R = Rlist, O = Olist, pen = pen.fac,
        sig.init = parameterr$sigma.fac[[j]],
        mont = mont, samples = samples, nb = n.b, no.it = no.it,
        s.p = s.p, delt = delt_sta,
        samp.train = samp.train, samp.test = samp.test, do.cv = TRUE)

      lambda2.1 <- matrix(gpu.fit$lambda, ncol = 1)
      sig       <- gpu.fit$sig
      deltaa    <- gpu.fit$deltaa
      delta.fac[j] <- deltaa

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
      hist.fac[[j]] <- list(deltaa = deltaa, pen.fac = pen.fac)
    }

    ################ M-step: structural (SEM) parameter estimation (GPU) #####
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
      if (length(which(aa == 1)) > 0) pen.sem <- pen.sem.list[[j]]

      ## ---- build SEM design g.eta / g.x (model-coupled, CPU) ----
      if (exists("et.ex")) rm(et.ex)
      if (exists("x.ex"))  rm(x.ex)
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
          g.eta <- lapply(1:mont, function(m) { lapply(1:samples, function(i1) 0) })
          et.ex <- "NULL"
        }
        if (length(obs) > 0) {
          if (any(co %in% model$var$observed)) {
            g.x <- list(); xx <- list()
            for (jj in 1:no.obs2) {
              fac <- which(model$mod$regression[[j]]$covariate == obs[jj])
              if (length(fac) == 0) next
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
            g.x <- lapply(1:samples, function(i1) 0); x.ex <- "NULL"
          }
        }
      }

      sig <- parameterr$sigma.sem[[lat.count]]
      if (!is.null(g.eta) | !is.null(g.x)) {
        ## assemble combined per-(m,sample) design g[[m]][[i1]] (with trimming)
        Dlist <- vector("list", mont); Rlist <- vector("list", mont); Olist <- vector("list", mont)
        for (m in 1:mont) {
          Dlist[[m]] <- vector("list", samples)
          Rlist[[m]] <- vector("list", samples)
          Olist[[m]] <- vector("list", samples)
          for (i1 in 1:samples) {
            geta <- if (is.null(g.eta)) NULL else g.eta[[m]][[i1]]
            gx   <- if (is.null(g.x))   NULL else g.x[[i1]]
            gmi <- cbind(geta, gx)
            if (exists("et.ex")) gmi <- gmi[, -1, drop = FALSE]
            if (exists("x.ex"))  gmi <- gmi[, -ncol(gmi), drop = FALSE]
            Dlist[[m]][[i1]] <- gmi
            Rlist[[m]][[i1]] <- as.vector(eta.l[[m]][[i1]][[lat.count]])
            Olist[[m]][[i1]] <- 0
          }
        }

        delt_sta <- c(0.001, 0.01, 0.05, 0.1) * log(samples)
        samp.train <- sort(sample(1:samples, round(0.7 * samples)))
        samp.test  <- sort((1:samples)[-samp.train])

        gpu.fit <- .fsem_gpu_penalized_cv_fit(
          D = Dlist, R = Rlist, O = Olist, pen = pen.sem,
          sig.init = parameterr$sigma.sem[[lat.count]],
          mont = mont, samples = samples, nb = n.b, no.it = no.it,
          s.p = s.p, delt = delt_sta,
          samp.train = samp.train, samp.test = samp.test, do.cv = TRUE)

        gamma2.1 <- matrix(gpu.fit$lambda, ncol = 1)
        sig      <- gpu.fit$sig
        deltaa   <- gpu.fit$deltaa
      } else {
        ## intercept-only: sigma update only (CPU, cheap)
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
        hist.sem[[lat.count]] <- list(deltaa = deltaa, pen.sem = pen.sem)
      } else {
        hist.sem[[lat.count]] <- list(deltaa = NA)
      }
      if (exists("et.ex")) rm(et.ex)
      if (exists("x.ex"))  rm(x.ex)
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

    param_prev <- params.eval
    sig_prev <- parameterr$sigma.error
    res.k = res.k + 1
    res.list[[res.k]] = list(parameters = parameters, parameters.eval = params.eval, timing = c(time.start = time.start, time.end = time.end))

    if (plot.progress) print(parameterr$sigma.error)

    if (!is.null(save.every)) {
      if (!is.null(save.every$steps) && (i %% save.every$steps) == 0) {
        dir_path = if (!is.null(save.every$dir)) save.every$dir else "."
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
