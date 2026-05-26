################################################################################
## run_MSE_tab1.R
##
## Stand-alone driver: open a fresh R session, set the working directory to
## the FSEM project root, and `source()` this file. It will compute the MSE
## tables from the saved Monte-Carlo results in outputs/ and print them.
##
## Usage (fresh R console):
##   setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
##   source("run_MSE_tab1.R")
##
## Optional fast check (averages only the first K runs per (N,M) cell):
##   Sys.setenv(MSE_MAX_PER_CELL = "4"); source("run_MSE_tab1.R")
##
## Returns (invisibly) a list with $fac and $sem; also assigns
## tabl.fac and tabl.sem in the global environment and writes
## outputs/MSE_tab1_fac.csv and outputs/MSE_tab1_sem.csv.
################################################################################

if (!file.exists("Core/MSE_tab1_fixed.R"))
  stop("Run this script from the FSEM project root (where 'Core/' lives).")

## Make sure we recompute from disk, not from a stale in-memory `results`.
if (exists("results", envir = .GlobalEnv)) rm("results", envir = .GlobalEnv)

source("Core/MSE_tab1_fixed.R")

cat("\nSaved: outputs/MSE_tab1_fac.csv, outputs/MSE_tab1_sem.csv\n")

invisible(list(fac = tabl.fac, sem = tabl.sem))
