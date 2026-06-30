## Aggregate _cov2_parts/*.part.rds -> Table 2 coverage by design x N x M.
setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
PARTDIR <- "_cov2_parts"
parts <- list.files(PARTDIR, pattern="\\.part\\.rds$", full.names=TRUE)
cat("parts found:", length(parts), "of 800\n")
res <- do.call(rbind, lapply(parts, readRDS))
# drop error rows (design NA)
nbad <- sum(is.na(res$design))
if(nbad>0) cat("WARNING:", nbad, "error/NA rows dropped\n")
res <- res[!is.na(res$design),]
saveRDS(res, "_cov_table2_raw.rds"); write.csv(res, "_cov_table2_raw.csv", row.names=FALSE)
agg <- aggregate(cbind(per1,per2,per3)~design+N+M, data=res, FUN=function(x) mean(x, na.rm=TRUE))
agg <- agg[order(agg$design, agg$N, agg$M),]
cat("\n============= TABLE 2 COVERAGE (lambda1/2/3) =============\n")
print(agg, row.names=FALSE, digits=4)
cat("\nrun counts per cell:\n")
print(table(res$design, res$N, res$M))
