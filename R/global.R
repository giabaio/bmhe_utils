#' This defines variables, mainly used with 'dplyr' and pipes that
#' throw a 'no visible binding for global variable' during
#' 'R CMD check' (before CRAN submission)
#'
#' @noRd
utils::globalVariables(c(
  ".",
  "n.used",
  "y",
  "Parameter",
  "num",
  "a",
  "b",
  "picker",
  "shape",
  "rate",
  "m",
  "required_packages",
  "variable",
  "value",
  "iteration",
  "chain"
)
)

