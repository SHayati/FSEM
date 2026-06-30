setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
FSEM_MSE_DEFINE_ONLY <- TRUE
suppressMessages({library(Matrix); library(tidyverse)})
source("Core/MSE_tab3.R")   # loads .mse_one_entry, .mse_aggregate, etc.
files <- list.files("outputs", pattern="^result_parameter_set1_regular\\.missing_[0-9]+_N_100_M_20\\.rds$", full.names=TRUE)
ent <- vector("list", length(files))
for(i in seq_along(files)){
  res <- readRDS(files[i])
  ent[[i]] <- .mse_one_entry(res, "model.fit", "estimation")
  rm(res); if(i%%20==0) gc()
}
tmpl <- readRDS(files[1])
tabs <- .mse_aggregate(ent, tmpl$model.fit, tmpl$estimation$result$params.estimated.eval)
cat("==== tabl.fac ====\n"); print(tabs$tabl.fac)
cat("==== tabl.sem ====\n"); print(tabs$tabl.sem)
saveRDS(tabs, "_s1_official.rds")
