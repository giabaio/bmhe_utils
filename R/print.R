#' Printing a jags object
#'
#' Printing a \code{jags} object
#'
#'
#' @param x an object of class `rjags', see \code{\link{jags}} for details
#' @param digits rounding for tabular output on the console (default is
#' to round to 1 decimal place)
#' @param intervals the quantiles for the posterior distribution to
#' be displayed in the summary statistics table
#' @param ... further arguments to \code{\link{print}}
#' @seealso \code{\link{jags}}
#' @author Gianluca Baio
#' @keywords print
print.rjags <- function(x, digits = 3, intervals = c(0.025, 0.25, 0.5, 0.75, 0.975), ...) {
  required_packages=c("R2jags")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("`", pkg, "` is required: install.packages('", pkg, "')")
    }
    if (requireNamespace(pkg, quietly = TRUE)) {
      if (!is.element(pkg, (.packages()))) {
        suppressMessages(suppressWarnings(attachNamespace(pkg)))
      }
    }
  }
  x <- x$BUGSoutput
  sims.matrix <- x$sims.matrix
  mu.vect <- apply(sims.matrix, 2, mean)
  sd.vect <- apply(sims.matrix, 2, sd)
  int.matrix <- apply(sims.matrix, 2, quantile, intervals)
  if (x$n.chains>1) {
    n.eff <- x$summary[, "n.eff"]
    Rhat <- x$summary[, "Rhat"]
  } else {
    n.eff <- Rhat <- NULL
  }
  summaryMatrix <- t(rbind(mu.vect, sd.vect, int.matrix, Rhat, n.eff))
  if(x$n.chains==1) {
    colnames(summaryMatrix) = c("mean","sd",paste0(intervals*100,"%"))
  } else {
    colnames(summaryMatrix) = c("mean","sd",paste0(intervals*100,"%"),"Rhat","n.eff")
  }

  rownameMatrix <- rownames(summaryMatrix)
  dev.idx <- match("deviance", rownameMatrix)
  if(any(!is.na(dev.idx))){
    summaryMatrix <- rbind(summaryMatrix[-dev.idx,], summaryMatrix[dev.idx,])
    rownames(summaryMatrix) <- c(rownameMatrix[-dev.idx], rownameMatrix[dev.idx])
  }

  if (!is.null(x$model.file))
    cat("Inference for Bugs model at \"", x$model.file, "\", ",
        sep = "")
  if (!is.null(x$program))
    cat("fit using ", x$program, ",", sep = "")
  cat("\n ", x$n.chains, " chains, each with ", x$n.iter, " iterations (first ",
      x$n.burnin, " discarded)", sep = "")
  if (x$n.thin > 1)
    cat(", n.thin =", x$n.thin)
  cat("\n n.sims =", x$n.sims, "iterations saved\n")
  print(round(summaryMatrix, digits), ...)
  if (x$n.chains > 1) {
    cat("\nFor each parameter, n.eff is a crude measure of effective sample size,")
    cat("\nand Rhat is the potential scale reduction factor (at convergence, Rhat=1).\n")
  }
  if (x$isDIC) {
    msgDICRule <- ifelse(x$DICbyR, "(using the rule, pD = var(deviance)/2)",
                         "(using the rule, pD = Dbar-Dhat)")
    cat(paste("\nDIC info ", msgDICRule, "\n", sep = ""))
    if (length(x$DIC) == 1) {
      cat("pD =", fround(x$pD, 1), "and DIC =", fround(x$DIC,
                                                       1))
    }
    else if (length(x$DIC) > 1) {
      print(round(x$DIC, 1))
    }
    cat("\nDIC is an estimate of expected predictive error (lower deviance is better).\n")
  }
  invisible(x)
}


#' Printing a bugs object
#'
#' Printing a \code{bugs} object
#'
#'
#' @param x an object of class `bugs', see \code{\link{bugs}} for details
#' @param digits.summary rounding for tabular output on the console (default is
#' to round to 1 decimal place)
#' @param intervals the quantiles for the posterior distribution to
#' be displayed in the summary statistics table
#' @param ... further arguments to \code{\link{print}}
#' @seealso \code{\link{bugs}}
#' @author Gianluca Baio
#' @keywords print
print.bugs <- function(x, digits = 3, intervals = c(0.025, 0.25, 0.5, 0.75, 0.975), ...) {
  required_packages=c("R2OpenBUGS")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("`", pkg, "` is required: install.packages('", pkg, "')")
    }
    if (requireNamespace(pkg, quietly = TRUE)) {
      if (!is.element(pkg, (.packages()))) {
        suppressMessages(suppressWarnings(attachNamespace(pkg)))
      }
    }
  }
  sims.matrix <- x$sims.matrix
  mu.vect <- apply(sims.matrix, 2, mean)
  sd.vect <- apply(sims.matrix, 2, sd)
  int.matrix <- apply(sims.matrix, 2, quantile, intervals)
  if (x$n.chains>1) {
    n.eff <- x$summary[, "n.eff"]
    Rhat <- x$summary[, "Rhat"]
  } else {
    n.eff <- Rhat <- NULL
  }
  summaryMatrix <- t(rbind(mu.vect, sd.vect, int.matrix, Rhat, n.eff))
  if(x$n.chains==1) {
    colnames(summaryMatrix) = c("mean","sd",paste0(intervals*100,"%"))
  } else {
    colnames(summaryMatrix) = c("mean","sd",paste0(intervals*100,"%"),"Rhat","n.eff")
  }

  rownameMatrix <- rownames(summaryMatrix)
  dev.idx <- match("deviance", rownameMatrix)
  if(any(!is.na(dev.idx))){
    summaryMatrix <- rbind(summaryMatrix[-dev.idx,], summaryMatrix[dev.idx,])
    rownames(summaryMatrix) <- c(rownameMatrix[-dev.idx], rownameMatrix[dev.idx])
  }

  if (!is.null(x$model.file))
    cat("Inference for Bugs model at \"", x$model.file, "\", ",
        sep = "")
  if (!is.null(x$program))
    cat("fit using ", x$program, ",", sep = "")
  cat("\n ", x$n.chains, " chains, each with ", x$n.iter, " iterations (first ",
      x$n.burnin, " discarded)", sep = "")
  if (x$n.thin > 1)
    cat(", n.thin =", x$n.thin)
  cat("\n n.sims =", x$n.sims, "iterations saved\n")
  print(round(summaryMatrix, digits), ...)
  if (x$n.chains > 1) {
    cat("\nFor each parameter, n.eff is a crude measure of effective sample size,")
    cat("\nand Rhat is the potential scale reduction factor (at convergence, Rhat=1).\n")
  }
  if (x$isDIC) {
    msgDICRule <- ifelse(x$DICbyR, "(using the rule, pD = var(deviance)/2)",
                         "(using the rule, pD = Dbar-Dhat)")
    cat(paste("\nDIC info ", msgDICRule, "\n", sep = ""))
    if (length(x$DIC) == 1) {
      cat("pD =", fround(x$pD, 1), "and DIC =", fround(x$DIC,
                                                       1))
    }
    else if (length(x$DIC) > 1) {
      print(round(x$DIC, 1))
    }
    cat("\nDIC is an estimate of expected predictive error (lower deviance is better).\n")
  }
  invisible(x)
  #print(x$BUGSoutput,...)
}

#' Rounding
#'
#' Rounds results
#'
#'
#' @param x an object of class `bugs', see \code{\link{bugs}} for details
#' @param digits number of digits to round
#' @seealso \code{\link{bugs}}
#' @author Gianluca Baio
#' @keywords print
#' @noRd
fround <- function(x, digits)
  format(round(x, digits), nsmall=digits)
