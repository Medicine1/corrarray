#' @title Correlation Arrays and 2-Sample Correlation Matrices
#'
#' @description This function creates a correlation array combining
#'   matrices for each level of a specified grouping variable.
#'   Given two levels of the grouping variable, this function creates a single
#'   correlation matrix displaying the individual triangular matrices on
#'   opposite sides of the principal diagonal.
#'
#' @usage corrarray(x, group = NULL, lower = NULL, upper = NULL,
#'   output = c("matrix", "array"),
#'   use = c("complete.obs", "everything", "all.obs",
#'   "na.or.complete", "pairwise.complete.obs"),
#'   method = c("pearson", "kendall", "spearman"))
#'
#' @param x a matrix or data frame. Variables can be quantitative or categorical.
#'
#' @param group the grouping variable name. If no group is specified (default),
#'   then a single correlation matrix for the entire sample will be generated.
#'
#' @param lower the level of the grouping variable to be placed in the
#'   lower triangular matrix. If no level is specified (default), then
#'   the first level in the data set is treated as the lower level.
#'
#' @param upper the level of the grouping variable to be placed in the
#'   upper triangular matrix. If no level is specified (default), then
#'   the second level in the data set is treated as the lower level.
#'
#' @param output If a group that has 2 or more levels is specified,
#'   then "\code{matrix}" returns the corresponding 2-sample correlation
#'   matrix (default), and "\code{array}" returns the correlation array.
#'
#' @param use an optional character string giving a method for
#'   computing correlations in the presence of missing values.
#'   This must be one of the strings "\code{complete.obs}" (default),
#'   "\code{everything}", "\code{all.obs}", "\code{na.or.complete}",
#'   or "\code{pairwise.complete.obs}".
#'   The default option removes rows with missing values from calculations.
#'
#' @param method a character string indicating which correlation coefficient
#'  is to be computed: "\code{pearson}" (default), "\code{kendall}",
#'  or "\code{spearman}".
#'
#' @details If multiple values are provided for \code{group},
#'   \code{lower}, or \code{upper}, then only the first value is used.
#'   Apart from the grouping variable, all other variables whose values are
#'   not numeric are removed from the correlation matrices. The grouping
#'   variable's values, even if numeric, are automatically treated as
#'   different levels.
#'
#' @return An array or matrix with numeric values from \code{-1} to
#'   \code{1} (inclusive) and with row and column names as the variable names.
#'
#' @export
#'
#' @seealso \code{\link{cor}} for further descriptions of the \code{use}
#'   and \code{method} parameters.
#'
#' @importFrom stats cor
#'
#' @examples
#' # All observations: 1-sample correlation matrix.
#' corrarray(iris)
#' # Stratify by the four species: 4-sample correlation array.
#' corrarray(iris, "Species", output = "array")
#' # Specify lower and upper samples: 2-sample correlation matrix.
#' corrarray(iris, "Species", lower = "setosa", upper = "virginica")
#'
corrarray <- function(x, group = NULL, lower = NULL, upper = NULL,
                      output = c("matrix", "array"),
                      use = c("complete.obs", "everything", "all.obs",
                              "na.or.complete", "pairwise.complete.obs"),
                      method = c("pearson", "kendall", "spearman")){

  # If group = NULL, return 1-sample correlation matrix.
  if (is.null(group)) {
    # Remove columns that are not numeric.
    new.x <- x[,unlist(lapply(x, FUN = is.numeric))]
    return(cor(new.x, use = "complete.obs", method = method))
  }

  # If multiple elements are specified, choose only first element.
  group <- group[1]
  lower <- lower[1]
  upper <- upper[1]

  # Identify whether group is a valid column name.
  if (!(group %in% names(x))) {
    stop(paste("Grouping variable '", group,
               "' is not found in data set.", sep = ""))
  }

  # If group is not referring to a factor variable, coerce it to factor.
  new.x <- x
  if (!is.factor(x[group])) {
    new.x[group] <- lapply(x[group] , as.factor)
  }

  # Remove columns that are not numeric vectors, keeping the grouping variable.
  new.x <- x[, (unlist(lapply(x, FUN = is.numeric)) | names(x) == group)]

  # Determine the number of groups and variables.
  grplevels <- levels(new.x[, group])
  ngroups <- length(grplevels)
  nvar <- ncol(new.x) - 1

  # Create a list of vectors with observation row numbers for each group.
  obs_groups <- list()
  for (i in 1:ngroups) {
    obs_groups[[i]] <- which(new.x[, group] == grplevels[i])
  }

  # Create a modified data frame with grouping variable column removed.
  new.x <- new.x[, names(new.x) != group]

  # Create a correlation array based on number of groups.
  corr.array <- array(dim = c(nvar, nvar, ngroups),
                      dimnames = list(names(new.x), names(new.x),
                                      Sample = grplevels))

  # Assign correlations to array for each group.
  for (i in 1:ngroups) {
    corr.array[,,i] <- cor(new.x[obs_groups[[i]],], use = "complete.obs",
                           method = method)
  }

  # Output array if output = "array".
  output <- match.arg(output)
  if (output == "array") {
    return(corr.array)
  }

  # Compute number of correlation combinations using formula 'ncor'='nvar'C2.
  ncor <- (factorial(nvar))/(2*factorial(nvar-2))

  # For 2-sample correlation matrix, define specified lower and upper groups.
  if (!is.null(lower) || !is.null(upper)) {
    if (!is.null(lower) && !(lower %in% grplevels)) {
      stop("Name of sample entered into 'lower' is not found in data set.")
    }

    else if (!is.null(upper) && !(upper %in% grplevels)) {
      stop("Name of sample entered into 'upper' is not found in data set.")
    }

    else if (!is.null(lower) && lower == upper) {
      stop("Names of samples entered into 'lower' and 'upper' are the same.")
    }

    else {
      lower <- which(grplevels == lower)
      upper <- which(grplevels == upper)
    }
  }

  else if (is.null(lower) && is.null(upper)) {
    lower <- 1
    upper <- 2
  }

  # Create matrix of proper dimensions. ####
  M <- matrix(data = NA, nrow = nvar, ncol = nvar,
              dimnames = list(Sample1=names(new.x),
                              Sample2=names(new.x)))

  # Comment with the factor level referred to by the sample numbers.
  comment(M) <- c(paste("Sample1 (lower triangular matrix) is '",
                        grplevels[lower], "' (n=",
                        sum(x[, group] == grplevels[lower]), ").",
                        sep = ""),
                  paste("Sample2 (upper triangular matrix) is '",
                        grplevels[upper], "' (n=",
                        sum(x[, group] == grplevels[upper]), ").",
                        sep = ""))

  # Fill in correlation matrix with corresponding correlations.
  for (i in 1:nvar) {
    for (j in 1:nvar) {
      if (j <= i) {
        M[i,j] <- corr.array[i,j,lower]
      }
      else if (j > i) {
        M[i,j] <- corr.array[i,j,upper]
      }
    }
  }

  # Output matrix and show the name of each sample.
  print(attr(x = M, which = "comment"))
  return(M)

}
