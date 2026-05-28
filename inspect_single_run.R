# Inspect single-run result from main_simulation_3.R
sr <- readRDS("single_run_results.rds")
cat("==== top-level names ====\n"); print(names(sr))
for (nm in names(sr)) {
  cat("\n----", nm, "----\n")
  if (is.list(sr[[nm]])) print(names(sr[[nm]])) else print(class(sr[[nm]]))
}
cat("\n==== sim.data structure ====\n")
str(sr$sim.data, max.level = 2)
cat("\n==== params.org.eval names ====\n")
print(names(sr$sim.data$params.org.eval))
cat("\n==== well-specified estimation names ====\n")
print(names(sr$params.estimated.well.specified))
print(names(sr$params.estimated.well.specified$result))
print(names(sr$params.estimated.well.specified$result$params.estimated.eval))
