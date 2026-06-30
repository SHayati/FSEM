setwd("C:/Users/shaya/Documents/VSCode/FSEM/FSEM")
r <- readRDS("outputs/result_parameter_set1_regular.missing_1_N_100_M_20.rds")
est <- r$estimation$result
cat("params.estimated names:", names(est$params.estimated), "\n")
cat("params.estimated.eval names:", names(est$params.estimated.eval), "\n\n")
# intercept: std and non-std, eval and non-eval
ie <- est$params.estimated.eval$intercept.std
inn <- est$params.estimated.eval$intercept      # non-std eval if exists
cat("intercept.std (eval)[[1]] mean:", round(mean(ie[[1]]),3), " head:", paste(round(head(ie[[1]],3),3),collapse=","),"\n")
if(!is.null(inn)) cat("intercept (eval)[[1]] mean:", round(mean(inn[[1]]),3),"\n")
op <- r$simulation$params.org.eval
cat("ORG intercept.std[[1]] mean:", round(mean(op$intercept.std[[1]]),3)," head:", paste(round(head(op$intercept.std[[1]],3),3),collapse=","),"\n")
# non-eval params.estimated intercept
cat("\nparams.estimated$intercept[[1]] (basis coefs) len:", length(est$params.estimated$intercept[[1]]),
    " vals:", paste(round(est$params.estimated$intercept[[1]],3),collapse=","),"\n")
# loadings non-std vs std
cat("\nlambda.param (non-std)[[2]] mean:", round(mean(est$params.estimated$lambda.param[[2]][[1]]),3),"\n")
cat("lambda.param.std[[2]] mean:", round(mean(est$params.estimated$lambda.param.std[[2]][[1]]),3),"\n")
# Does coef.fac.std (eval) differ from lambda.param.std (eval)?
cat("coef.fac.std.eval[[2]][[1]] mean:", round(mean(est$params.estimated.eval$coef.fac.std[[2]][[1]]),3),"\n")
# org loadings non-eval? Check the org raw lambda
cat("\nORG coef.fac.std.eval[[2]][[1]] mean:", round(mean(op$coef.fac.std[[2]][[1]]),3),"\n")
# sigma error
cat("\nsigma.error est:", round(unlist(est$params.estimated$totparam$sigma.error),3),"\n")
cat("sigma.error org:", round(unlist(op$sigma.error),3),"\n")
