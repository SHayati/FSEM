## Verify the gamma-stuck-at-zero hypothesis on the saved single run
sr <- readRDS("single_run_results.rds")

cat("=== TRUE gamma(t) eta1<-eta2  (head/tail/range) ===\n")
gt <- sr$sim.data$params.org.eval$coef.sem.std[["eta1"]][["eta2"]]
cat("length:", length(gt), "\n")
cat("range: [", round(min(gt), 4), ",", round(max(gt), 4), "]\n")
cat("mean(|true|): ", round(mean(abs(gt)), 4), "\n")

cat("\n=== ESTIMATED gamma(t) (well-specified) ===\n")
ge <- sr$params.estimated.well.specified$result$params.estimated.eval$coef.sem.std[["eta1"]][["eta2"]]
cat("length:", length(ge), "\n")
cat("range: [", signif(min(ge), 4), ",", signif(max(ge), 4), "]\n")
cat("mean(|est|): ", signif(mean(abs(ge)), 4), "\n")

cat("\n=== RAW gamma.param vector (before evaluation, well-specified) ===\n")
gp <- sr$params.estimated.well.specified$result$params.estimated$totparam$gamma.param
for (i in seq_along(gp)) {
  cat(sprintf("  gamma.param[[%d]] (len=%d): max|.|= %.3e   nonzero entries: %d\n",
              i, length(gp[[i]]), max(abs(gp[[i]])), sum(abs(gp[[i]]) > 1e-10)))
  if (length(gp[[i]]) <= 20) print(round(gp[[i]], 6))
}

cat("\n=== Trace tolerance (convergence) ===\n")
tol <- sr$params.estimated.well.specified$tolerance
cat("steps =", sr$params.estimated.well.specified$steps, "\n")
cat("tol (last 10): ", round(tail(tol$tol, 10), 5), "\n")
cat("tol_sig (last 10): ", round(tail(tol$tol_sig, 10), 5), "\n")
