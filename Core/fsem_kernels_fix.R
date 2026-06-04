## Core/fsem_kernels_fix.R
## ---------------------------------------------------------------------------
## Robust, collision-free loader for the optional Rcpp kernels used by
## em.estimation.optimised.
##
## Why this exists:
##   The default loader in FSEM_optimised.R calls Rcpp::sourceCpp() lazily on
##   the first EM iteration, using Rcpp's *shared* default cache directory.
##   When several R processes build the same kernel at the same time -- e.g.
##   the parallel workers in main_simulation_2_parallel.R, or two scripts run
##   back-to-back while a DLL is still mapped -- they race to write the same
##   temporary .dll. One process locks the file, the linker step in the others
##   fails, and R reports:
##       "Rcpp kernels unavailable (Error 322 occurred building shared
##        library.); using pure-R fallback."
##   The results are still correct (the pure-R path is numerically identical),
##   just slower.
##
## What this fixes it:
##   Compile the kernels ONCE per process into a private, PID-specific cache
##   directory, so no two processes ever touch the same build artefacts. Then
##   set the package flags so the in-pipeline loader short-circuits and uses
##   these functions directly.
##
## Usage: source this AFTER Core/FSEM_optimised.R, e.g.
##   source("Core/FSEM.R")
##   source("Core/simualtion_utilities.R")
##   source("Core/FSEM_optimised.R")
##   source("Core/fsem_kernels_fix.R")
## In parallel scripts, source it inside each worker (after FSEM_optimised.R).
##
## This file does NOT modify Core/FSEM_optimised.R or Core/FSEM.R.
## ---------------------------------------------------------------------------

.fsem_load_kernels_robust <- function(verbose = FALSE) {
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    assign(".fsem_have_rcpp", FALSE, envir = globalenv())
    return(invisible(FALSE))
  }
  candidates <- c(
    file.path("Core", "fsem_kernels.cpp"),
    "fsem_kernels.cpp",
    file.path(getwd(), "Core", "fsem_kernels.cpp")
  )
  cpp <- candidates[file.exists(candidates)][1]
  if (is.na(cpp)) {
    assign(".fsem_have_rcpp", FALSE, envir = globalenv())
    return(invisible(FALSE))
  }

  ## ensure the target env exists (normally created by FSEM_optimised.R)
  if (!exists(".fsem_kernel_env", envir = globalenv(), inherits = FALSE)) {
    assign(".fsem_kernel_env", new.env(parent = globalenv()), envir = globalenv())
  }
  kenv <- get(".fsem_kernel_env", envir = globalenv())

  ## private per-process cache dir -> parallel workers never collide
  cdir <- file.path(tempdir(), sprintf("fsem_kernels_%d", Sys.getpid()))
  dir.create(cdir, showWarnings = FALSE, recursive = TRUE)

  ok <- tryCatch({
    Rcpp::sourceCpp(cpp, env = kenv, cacheDir = cdir, rebuild = TRUE,
                    verbose = verbose)
    TRUE
  }, error = function(e) {
    message("fsem_kernels_fix: kernel precompile failed (",
            conditionMessage(e), "); using pure-R fallback.")
    FALSE
  })

  assign(".fsem_have_rcpp", ok, envir = globalenv())
  if (isTRUE(ok) && verbose)
    message("fsem_kernels_fix: Rcpp kernels compiled and loaded (pid ",
            Sys.getpid(), ").")
  invisible(ok)
}

## run once at source time
.fsem_load_kernels_robust()
