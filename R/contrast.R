contrast <- function(X, abar, nodes, contrastType = "difference",
                     measure = mean, cond = NULL, cilevel = 0.95, ...) {
  # --------------------------- Check input arguments --------------------------
  if ("simulated.data" %in% names(X)) {
    setup <- "gformula"
    if (is.null(X[["simulated.data"]])) {
      stop("Set 'ret=TRUE' in gformula().")
    }
  } else {
    setup <- "sgf"
  }

  # Survival setting: transform matrices to dataframes
  if (setup == "gformula" && is.matrix(X$simulated.data[[1]][[1]])) {
    X$simulated.data <- lapply(X$simulated.data, function(blevel) {
      lapply(blevel, function(m) as.data.frame(m))
    })
  }

  if (setup == "gformula" &&
        !all(nodes %in% names(X$simulated.data[[1]][[1]]))) {
    stop("'nodes' contains names that are not present in the data.")
  }

  if (is.character(contrastType) && !(length(nodes) %in% c(1, 2))) {
    stop("'nodes' must be of length 1 or 2.")
  }

  if (!(
    (is.character(contrastType) && contrastType %in%
       c("difference", "ratio", "oddsratio")) ||
      (is.function(contrastType))
  )) {
    stop("'contrastType' must be either one of the character strings:\n",
         "'difference', 'ratio', 'oddsratio', or a function.")
  }

  if (!(is.numeric(cilevel) && length(cilevel) == 1 &&
          cilevel > 0 && cilevel < 1)) {
    stop("'cilevel' must be a single numeric value between 0 and 1.")
  }

  if (setup == "sgf") {
    if (!identical(measure, base::mean)) {
      warning(
        "When calculating a contrast based on sgf(), 'measure' must be 'mean'.",
        " Resetting 'measure' to 'mean'."
      )
      measure <- base::mean
    }
    if (!is.null(cond)) {
      warning(
        "When calculating a contrast based on sgf(), conditioning with 'cond' ",
        "is not possible. Resetting 'cond' to NULL."
      )
      cond <- NULL
    }
  }

  if (!is.function(measure)) {
    stop("'measure' must be a function.")
  }

  if (!is.null(cond)) {
    if (is.language(cond)) {
      condVars <- all.vars(cond)
      numCond <- 1
    } else if (is.list(cond) && all(vapply(cond, is.language, logical(1)))) {
      condVars <- unique(unlist(lapply(cond, all.vars)))
      numCond <- length(cond)
    } else {
      stop(
        "'cond' must be a quoted expression or a list of quoted expressions,
         e.g., quote(sex == 1) or list(quote(sex == 1), quote(sex == 0))."
      )
    }
    condNA <- setdiff(condVars, names(X$simulated.data[[1]][[1]]))
    if (length(condNA) > 0) {
      stop(sprintf(
        "The following variables used in 'cond' are not present in the data:
        %s.\n", paste(condNA, collapse = ", ")
      ))
    }
    if (is.character(contrastType) && !(numCond %in% c(1, 2))) {
      stop("'cond' must be of length 1 or 2 for 'contrastType = %s'.",
           contrastType)
    }
  }

  # -------------------- Handle different intervention types -------------------
  if (is.matrix(abar)) {
    abarStr <- apply(abar, 1, paste, collapse = ",")
    setupStr <- apply(X$setup$abar, 1, paste, collapse = ",")
    aidx <- match(abarStr, setupStr)
    if (anyNA(aidx)) {
      stop(sprintf(
        "'abar' contains rows that were not evaluated in %s().", setup
      ))
    }
  } else if (is.vector(abar) && is.numeric(abar)) {
    if (X$setup$i.type != "individual") {
      aidx <- match(abar, X$setup$abar)
    } else {
      abarStr <- as.character(abar)
      setupStr <- sapply(X$setup$abar, function(val) unique(colnames(val)))
      aidx <- match(abarStr, setupStr)
    }
    if (anyNA(aidx)) {
      stop(sprintf(
        "'abar' contains values that were not evaluated in %s().", setup
      ))
    }
  } else if ("natural" %in% abar) {
    if (!("natural" %in% X$setup$abar)) {
      stop("'gformula()' was not evaluated with 'abar' = 'natural'.")
    }
  } else if (is.list(abar)) {
    abarStr <- vapply(
      abar, function(mat) paste(sort(colnames(mat)), collapse = ","),
      FUN.VALUE = character(1)
    )
    setupStr <- vapply(
      X$setup$abar, function(mat) paste(sort(colnames(mat)), collapse = ","),
      FUN.VALUE = character(1)
    )
    aidx <- match(abarStr, setupStr)
    if (anyNA(aidx)) {
      stop(sprintf(
        "'abar' contains elements that were not evaluated in %s().", setup
      ))
    }
  } else {
    stop("'abar' incorrectly specified.")
  }

  if (is.character(contrastType)) {
    if (is.null(cond) || numCond == 1) {
      if (!("natural" %in% abar) && (!(length(aidx) %in% c(1, 2)) ||
                                       length(nodes) + length(aidx) != 3)) {
        stop(
          "'nodes' and 'abar' must be of compatible lengths:\n",
          "- 'nodes' can be length 1 only if 'abar' contains 2 interventions\n",
          "- 'abar' can contain 1 intervention only if 'nodes' is of length 2"
        )
      }
    } else {
      if (
        !("natural" %in% abar) &&
          (!length(aidx) %in% c(1, 2) ||
             length(nodes) + length(aidx) + numCond != 4)
      ) {
        stop(
          "'nodes', 'abar', and 'cond' must be of compatible lengths:\n",
          "- 'nodes' can be length 1 if 'abar' contains 2 interventions\n",
          "- 'abar' can contain 1 intervention if 'nodes' is of length 2\n",
          "- 'nodes' and 'abar' can both be length 1 if 'cond' is of length 2"
        )
      }
    }
  } else {
    if (!("natural" %in% abar)) {
      expectedArgs <- if (is.null(cond)) {
        length(nodes) * length(aidx)
      } else {
        length(nodes) * length(aidx) * numCond
      }
      actualArgs <- length(setdiff(names(formals(contrastType)), "..."))

      if (expectedArgs != actualArgs) {
        stop(sprintf(
          "`contrastType` must be a function of %d arguments (as defined by
          'nodes', 'abar', and 'cond'), but it currently has %d.",
          expectedArgs, actualArgs
        ))
      }
    }
  }

  if (setup == "sgf") {
    tidx <- match(nodes, X$setup$Ynodes)
    if (anyNA(tidx)) {
      stop("'nodes' must only contain Ynodes that were evaluated in sgf().")
    }
  }

  # ----------------------------- Evaluate contrast ----------------------------
  boots <- NULL
  B <- if (setup == "sgf") ncol(X$boot.results) else length(X$simulated.data)
  for (b in 1:B) {
    if (setup == "gformula") {
      out <- NULL
      # Extract the column measures corresponding to the counterfactuals defined
      # by 'abar', 'nodes', and 'cond':
      if (!("natural" %in% abar)) {
        for (node in nodes) {
          for (a in aidx) {
            df <- X$simulated.data[[b]][[a]]
            if (!is.null(cond)) {
              if (numCond == 1) {
                df_cond <- subset(df, eval(cond, envir = df,
                                           enclos = parent.frame()))
                out <- c(out, measure(df_cond[[node]], ...))
              } else {
                for (cnd in cond) {
                  df_cond <- subset(df, eval(cnd, envir = df,
                                             enclos = parent.frame()))
                  out <- c(out, measure(df_cond[[node]], ...))
                }
              }
            } else {
              out <- c(out, measure(df[[node]], ...))
            }
          }
        }
      } else {
        if (length(nodes) > 1 || (!is.null(cond) && numCond > 1)) {
          for (node in nodes) {
            df <- X$simulated.data[[b]][[1]]
            if (!is.null(cond)) {
              if (numCond == 1) {
                df_cond <- subset(df, eval(cond, envir = df,
                                           enclos = parent.frame()))
                out <- c(out, measure(df_cond[[node]], ...))
              } else {
                for (cnd in cond) {
                  df_cond <- subset(df, eval(cnd, envir = df,
                                             enclos = parent.frame()))
                  out <- c(out, measure(df_cond[[node]], ...))
                }
              }
            } else {
              out <- c(out, measure(df[[node]], ...))
            }
          }
        } else {
          dfnatural <- X$simulated.data[[b]][[1]]
          dfobserved <- X$observed.data[[b]]
          if (!is.null(cond)) {
            dfnatural <- subset(
              dfnatural, eval(cond, envir = dfnatural, enclos = parent.frame())
            )
            dfobserved <- subset(
              dfobserved, eval(cond, envir = dfobserved,
                               enclos = parent.frame())
            )
          }
          out[1] <- measure(dfnatural[[nodes]], ...)
          out[2] <- measure(dfobserved[[nodes]], ...)
        }
      }
    } else {
      if (is.matrix(abar)) {
        abarSetupLength <- nrow(X$setup$abar)
      } else {
        abarSetupLength <- length(X$setup$abar)
      }
      out <- X$boot.results[(tidx - 1) * abarSetupLength + aidx, b]
    }
    if (is.character(contrastType)) {
      if (contrastType == "difference") {
        res <- out[1] - out[2]
      } else if (contrastType == "ratio") {
        res <- out[1] / out[2]
      } else if (contrastType == "oddsratio") {
        res <- (out[1] * (1 - out[2])) / (out[2] * (1 - out[1]))
      }
    } else {
      res <- do.call(contrastType, as.list(out))
    }
    boots <- c(boots, res)
    if (b == 1) {
      if (setup == "sgf") {
        out <- X$results[(tidx - 1) * abarSetupLength + aidx, "psi"]
        if (is.character(contrastType)) {
          if (contrastType == "difference") {
            res <- out[1] - out[2]
          } else if (contrastType == "ratio") {
            res <- out[1] / out[2]
          } else if (contrastType == "oddsratio") {
            res <- (out[1] * (1 - out[2])) / (out[2] * (1 - out[1]))
          }
        } else {
          res <- do.call(contrastType, as.list(out))
        }
      }
      pointEst <- list(counterfactuals = out,
                       contrast = res)
    }
  }

  # ------------------------------- Bootstrap CI -------------------------------
  fullboots <- if (setup == "gformula") boots[-1] else boots
  cibounds <- quantile(fullboots, probs = c((1 - cilevel) / 2,
                                            1 - (1 - cilevel) / 2))
  varhat <- var(fullboots)

  # ------------------------------- Return value -------------------------------
  returnObj <- list("counterfactuals" = pointEst[["counterfactuals"]],
                    "contrast" = pointEst[["contrast"]],
                    "ciContrast" = cibounds,
                    "B" = length(fullboots),
                    "varContrast" = varhat)

  numNodes <- length(nodes)
  numCondAlt <- if (!is.null(cond)) numCond else 1
  numAbar <- length(returnObj[["counterfactuals"]]) / (numNodes * numCondAlt)
  if (sum(c(numNodes, numCondAlt, numAbar) > 1) == 1) {
    if (numNodes > 1) {
      names(returnObj[["counterfactuals"]]) <- nodes
    } else if (numCondAlt > 1) {
      names(returnObj[["counterfactuals"]]) <- cond
    } else {
      if ("natural" %in% abar) {
        names(returnObj[["counterfactuals"]]) <- c("natural", "observed")
      } else if (is.matrix(abar)) {
        names(returnObj[["counterfactuals"]]) <- paste(
          "a =", apply(abar, 1, function(x) paste(round(x, 4), collapse = ","))
        )
      } else if (is.list(abar)) {
        names(returnObj[["counterfactuals"]]) <- paste("a =", vapply(
          abar, function(mat) colnames(mat)[1], FUN.VALUE = character(1)
        ))
      } else {
        names(returnObj[["counterfactuals"]]) <- paste("a =", round(abar, 4))
      }
    }
  } else {
    if ("natural" %in% abar) {
      abarNames <- "natural"
    } else if (is.matrix(abar)) {
      abarNames <- paste(
        "a =", apply(abar, 1, function(x) paste(round(x, 4), collapse = ","))
      )
    } else if (is.list(abar)) {
      abarNames <- paste("a =", vapply(
        abar, function(mat) colnames(mat)[1], FUN.VALUE = character(1)
      ))
    } else {
      abarNames <- paste("a =", round(abar, 4))
    }

    combNames <- NULL
    for (nodeName in nodes) {
      for (abarName in abarNames) {
        if (numCondAlt == 1) {
          combNames <- c(combNames, paste(nodeName, abarName, sep = ", "))
        } else {
          for (condName in cond) {
            depCondName <- paste(deparse(condName), collapse = "")
            combNames <- c(
              combNames, paste(nodeName, abarName, depCondName, sep = ", ")
            )
          }
        }
      }
    }
    names(returnObj[["counterfactuals"]]) <- combNames
  }
  if (is.character(contrastType)) {
    names(returnObj[["contrast"]]) <- contrastType
  } else {
    tryCatch(
      {
        nm <- paste(deparse(substitute(contrastType)), collapse = "")
        names(returnObj[["contrast"]]) <- nm
      },
      error = function(e) {
        names(returnObj[["contrast"]]) <- "customcontrast"
      }
    )
  }

  class(returnObj) <- "contrastResult"
  return(returnObj)
}