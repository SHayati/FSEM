## Independent worker for Table 2 coverage. Args: worker_id n_workers
## Processes its stripe of files, writes part rds, skips done. Crash-isolated.
setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
args <- commandArgs(trailingOnly=TRUE)
wid <- as.integer(args[1]); nw <- as.integer(args[2])
source("_cov_fun.R")
if(!dir.exists(PARTDIR)) dir.create(PARTDIR)

files <- list.files("outputs", pattern="^result_parameter_set1_(regular|irregular)_[0-9]+_N_[0-9]+_M_[0-9]+\\.rds$", full.names=TRUE)
mine <- files[((seq_along(files)-1) %% nw) + 1 == wid]
logf <- sprintf("_cov2_w%02d.log", wid)
cat(sprintf("worker %d/%d: %d files\n", wid, nw, length(mine)), file=logf)

for(f in mine){
  partf <- file.path(PARTDIR, paste0(basename(f), ".part.rds"))
  if(file.exists(partf)) next
  set.seed(1000 + which(files==f))
  out <- tryCatch(cov_one_file(f), error=function(e){
    d<-data.frame(design=NA,N=NA,M=NA,per1=NA,per2=NA,per3=NA,file=basename(f),stringsAsFactors=FALSE)
    attr(d,"err")<-conditionMessage(e); d })
  saveRDS(out, partf)
  cat(sprintf("[%s] done %s\n", format(Sys.time(),"%H:%M:%S"), basename(f)), file=logf, append=TRUE)
  rm(out); gc()
}
cat("worker finished\n", file=logf, append=TRUE)
