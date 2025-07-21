print.contrastResult <- function(x, ...) {
  cf <- x$counterfactuals
  cType <- names(x$contrast)
  est <- x$contrast
  ci <- x$ciContrast
  B <- x$B

  cat("------------------------------------------------------------\n")
  cat("The estimated counterfactual outcomes are:\n")
  for (i in seq_along(cf)) {
    cat(sprintf("    %-12s %12.4f\n", names(cf)[i], cf[i]))
  }
  cat("\n")

  if (length(cf) == 2) {
    cat(sprintf("The estimated %s between %s and %s is:\n",
                cType, names(cf)[1], names(cf)[2]))
  } else {
    cat(sprintf("The estimated %s is:\n", cType))
  }
  cat(sprintf("    %12.4f\n", est))
  cat("\n")

  if (B > 0) {
    cat(sprintf("Based on B = %d successful bootstrap samples, the\n", B))
    cat(sprintf("estimated confidence interval for the %s is:\n", cType))
    cat(sprintf("    %-8s %16s\n", names(ci)[1], names(ci)[2]))
    cat(sprintf("    %8.4f %16.4f\n", ci[1], ci[2]))
  } else {
    cat("To estimate a confidence interval, specify the number\n")
    cat("of bootstrap samples B > 1 in 'gformula'.\n")
  }
  cat("------------------------------------------------------------\n")

  invisible(x)
}