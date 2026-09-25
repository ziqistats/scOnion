#' The Onion object
#'
#' An `Onion` is the fitted scOnion model returned by [Onion()]. It is a plain
#' named list carrying an S3 class attribute, so every element remains reachable
#' with `$`. The class adds `print` and `summary` methods, lets the other
#' functions check that they were handed a real fit, and keeps the element
#' names used before the package was renamed working (see the section on legacy
#' names).
#'
#' The elements are named after the symbols of the scOnion Methods. The picture
#' behind the name: the two denoised modalities are the two halves of an onion,
#' each a nested set of rings, one ring per singular component. The SVD of the
#' small core where the halves meet gives the paired coordinates `U` and `V`,
#' and laying the two halves side by side gives the cross-section `W`.
#'
#' @section Elements:
#' \describe{
#'   \item{`W`}{`n` x `(r3 + r4)` matrix. The joint cell embedding
#'     \eqn{W = [\hat{U} | \hat{V}]}. Rows are cells, named by barcode; columns
#'     are named `onion_1`, `onion_2`, ...}
#'   \item{`U`, `V`}{`n` x `r3` and `n` x `r4` matrices. The paired cell
#'     embeddings, one per modality. Their columns are matched one to one
#'     through the SVD, which is what makes the per-cell concordance score
#'     meaningful.}
#'   \item{`stdev.cell`}{singular values of the cell embedding.}
#'   \item{`P`, `Q`}{`g` x `r5` and `p` x `r5` matrices. Genome-wide gene and
#'     second-modality feature loadings.}
#'   \item{`stdev.feature`}{singular-value weights \eqn{d} of the feature
#'     embedding. These enter the linkage score.}
#'   \item{`stdev.feature.full`}{the full feature spectrum, kept so that the
#'     denoiser figure can be redrawn from a saved fit; see [onionSpectrum()].}
#'   \item{`per.chr`}{named list, one entry per chromosome, each
#'     `list(P, Q, stdev)`. This is the per-chromosome decomposition that
#'     [OnionLinkageTable()] scores, so every linkage is within-chromosome.}
#'   \item{`X.svd`, `Y.svd`}{the truncated SVDs of the two denoised modalities.
#'     Retained because [OnionLinkageTable()] reconstructs the denoised feature
#'     vectors from them to compute the jackknife standard error.}
#'   \item{`input.param`}{the thresholds supplied, the ranks actually selected,
#'     and the versions of the software that produced the fit.}
#' }
#'
#' Optional elements are added in place by downstream functions:
#' `concordance` and `simscore.*` by [OnionConcordance()], `UMAPs` by
#' [OnionUMAP()], and `NN.graph` by [OnionFindNeighbors()].
#'
#' @section Legacy names:
#' Until version 0.99.1 the package was called POLARIS and the fit was a
#' `skymap` with elements `skymap.cell`, `skymap.U`, `skymap.V`, `skymap.P`,
#' `skymap.Q` and `skymap.feature.chr`. Fits saved under those names, including
#' fits made by the pre-package script, still work: every function upgrades
#' them on entry, and `fit$skymap.cell` on an `Onion` returns `W` with a
#' warning, shown once per session and name. Convert a saved fit once with
#' [asOnion()] to stop the warnings. The legacy names will be removed in
#' version 1.0.0.
#'
#' @name Onion-class
#' @aliases Onion-object
#' @seealso [asOnion()], [scOnion-deprecated]
NULL

## Element names used before the package was renamed, and their current names.
## Read by asOnion(), by the `$` and `[[` methods, and by .resolve_slots().
.legacy_slots <- c(skymap.cell = "W", skymap.U = "U", skymap.V = "V",
                   skymap.P = "P", skymap.Q = "Q",
                   skymap.feature.chr = "per.chr")

## Once-per-session bookkeeping for the legacy-name warnings and message.
.legacy_env <- new.env(parent = emptyenv())

#' @noRd
.warn_legacy <- function(old, new) {
  key <- paste0("slot:", old)
  if (isTRUE(.legacy_env[[key]])) return(invisible(NULL))
  assign(key, TRUE, envir = .legacy_env)
  warning(sprintf(paste0("`%s` is a POLARIS-era name; use `%s`. Legacy names ",
                         "still work in scOnion 0.99.x and will be removed in ",
                         "1.0.0. Run asOnion() on a saved fit to update it. ",
                         "This warning is shown once per session."),
                  old, new), call. = FALSE)
  invisible(NULL)
}

#' @noRd
.new_Onion <- function(x) {
  structure(x, class = c("Onion", "list"))
}

#' Is `x` a fit in the old layout, or an unclassed fit?
#' @noRd
.is_legacy_fit <- function(x) {
  is.list(x) && !inherits(x, "Onion") &&
    (inherits(x, "skymap") ||
       !is.null(.subset2(x, "skymap.cell")) || !is.null(.subset2(x, "W")))
}

#' Convert a saved POLARIS fit to the Onion layout
#'
#' Renames the elements of a fit made by POLARIS 0.99.0 or 0.99.1, or by the
#' pre-package script, to the names an [Onion][Onion-class] uses, and sets the
#' class. No matrix is copied, only the names change. An `Onion` is returned
#' unchanged, so the call is safe to apply everywhere.
#'
#' @param x A fit: an `Onion`, a legacy `skymap`, or the bare list the
#'   pre-package script returned.
#'
#' @return An `Onion`.
#'
#' @details
#' The renames are `skymap.cell` to `W`, `skymap.U` to `U`, `skymap.V` to `V`,
#' `skymap.P` to `P`, `skymap.Q` to `Q` and `skymap.feature.chr` to `per.chr`,
#' together with the derived keys under `UMAPs` (`<slot>.UMAP`) and `NN.graph`
#' (`<slot>.nb`). `input.param$polaris.version` is copied to `onion.version`
#' and left in place, since it records the software that made the fit.
#'
#' Every function in the package applies this conversion on entry when handed a
#' legacy fit, so calling it yourself is only needed to silence the message and
#' to resave the fit in the current layout.
#'
#' @examples
#' set.seed(1)
#' sim <- onionSimulate(n = 150, g = 50, p = 60, r = 4)
#' fit <- Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE)
#' old <- unclass(fit)
#' names(old)[names(old) == "W"] <- "skymap.cell"   # as a saved 0.99.1 fit looks
#' identical(asOnion(old)$W, fit$W)
#'
#' @seealso [Onion-class], [scOnion-deprecated]
#' @export
asOnion <- function(x) {
  if (inherits(x, "Onion")) return(x)
  if (!.is_legacy_fit(x))
    stop(paste0("`x` must be an Onion fit returned by Onion(), or a fit saved ",
                "by POLARIS 0.99.x or the pre-package script."), call. = FALSE)
  nm <- names(x)
  hit <- nm %in% names(.legacy_slots)
  nm[hit] <- unname(.legacy_slots[nm[hit]])
  names(x) <- nm
  for (sl in c("UMAPs", "NN.graph")) {
    l <- .subset2(x, sl)
    if (is.null(l) || is.null(names(l))) next
    kn <- names(l)
    for (old in names(.legacy_slots)) {
      pat <- paste0("^", gsub(".", "\\.", old, fixed = TRUE), "(\\.UMAP|\\.nb)$")
      kn <- sub(pat, paste0(.legacy_slots[[old]], "\\1"), kn)
    }
    names(l) <- kn
    x[[sl]] <- l
  }
  ip <- .subset2(x, "input.param")
  if (is.list(ip) && !is.null(ip$polaris.version) && is.null(ip$onion.version)) {
    ip$onion.version <- ip$polaris.version
    x[["input.param"]] <- ip
  }
  .new_Onion(x)
}

#' Check, and if needed upgrade, the fit handed to a verb
#'
#' Returns the fit, converted to the current layout when it was saved under the
#' POLARIS names, so callers write `object <- .check_Onion(object, ...)`.
#' @noRd
.check_Onion <- function(object, need = character()) {
  if (!inherits(object, "Onion")) {
    if (!.is_legacy_fit(object))
      stop("`object` must be an Onion fit returned by Onion().", call. = FALSE)
    if (!isTRUE(.legacy_env[["upgraded"]])) {
      assign("upgraded", TRUE, envir = .legacy_env)
      message(paste0("Upgrading a POLARIS-era fit to the Onion layout. Run ",
                     "asOnion() on saved fits to do this once. ",
                     "(Shown once per session.)"))
    }
    object <- asOnion(object)
  }
  missing <- need[vapply(need, function(n) is.null(.subset2(object, n)),
                         logical(1))]
  if (length(missing))
    stop(sprintf(paste0("This fit has no %s. Re-run Onion(), or the function ",
                        "that produces it."),
                 paste(missing, collapse = ", ")), call. = FALSE)
  object
}

## `$` and `[[` on an Onion accept the POLARIS-era element names, so that a
## script written against `fit$skymap.cell` keeps working (with a warning)
## instead of silently receiving NULL. .subset2() is the internal default, so
## there is no recursion; exact = FALSE keeps the partial matching that `$` on a
## plain list has.

#' @export
`$.Onion` <- function(x, name) {
  new <- .legacy_slots[name]
  if (!is.na(new)) {
    .warn_legacy(name, new)
    name <- unname(new)
  }
  .subset2(x, name, exact = FALSE)
}

#' @export
`[[.Onion` <- function(x, i, ...) {
  if (is.character(i) && length(i) == 1L) {
    new <- .legacy_slots[i]
    if (!is.na(new)) {
      .warn_legacy(i, new)
      i <- unname(new)
    }
  }
  .subset2(x, i, ...)
}

#' @method print Onion
#' @export
print.Onion <- function(x, ...) {
  n <- nrow(x$W)
  cat(sprintf("An Onion: scOnion fit on %s cells\n", format(n, big.mark = ",")))
  cat(sprintf("  cross-section W  %s x %d   (halves U: %d, V: %d)\n",
              format(n, big.mark = ","), ncol(x$W),
              ncol(x$U), ncol(x$V)))
  if (!is.null(x$P))
    cat(sprintf("  loadings         P  %s x %d   Q  %s x %d\n",
                format(nrow(x$P), big.mark = ","), ncol(x$P),
                format(nrow(x$Q), big.mark = ","), ncol(x$Q)))
  if (!is.null(x$per.chr)) {
    ok <- sum(!vapply(x$per.chr, is.null, logical(1)))
    cat(sprintf("  per-chromosome   %d of %d chromosomes decomposed\n",
                ok, length(x$per.chr)))
  }
  ip <- x$input.param
  if (!is.null(ip) && !is.null(ip$T1)) {
    conv <- if (is.null(ip$reconcile.ranks)) ""
            else if (isTRUE(ip$reconcile.ranks)) "   (ranks reconciled)"
            else "   (ranks kept per modality)"
    cat(sprintf("  thresholds       T1=%g T2=%g T3=%g%s\n",
                ip$T1, ip$T2, ip$T3, conv))
  }
  extra <- setdiff(names(x), c("W", "U", "V", "stdev.cell", "P", "Q",
                               "stdev.feature", "stdev.feature.full", "per.chr",
                               "X.svd", "Y.svd", "input.param"))
  if (length(extra))
    cat(sprintf("  also present     %s\n", paste(extra, collapse = ", ")))
  invisible(x)
}

#' @method summary Onion
#' @export
summary.Onion <- function(object, ...) {
  print(object)
  if (!is.null(object$concordance)) {
    q <- stats::quantile(object$concordance, c(0, .25, .5, .75, 1), na.rm = TRUE)
    cat("\nconcordance (rescaled to [0,1]):\n")
    print(round(q, 4))
  }
  invisible(object)
}

## A fit saved by POLARIS 0.99.x carries class "skymap". Printing it upgrades
## it on the fly, so a saved object still displays without user action.

#' @method print skymap
#' @export
print.skymap <- function(x, ...) {
  print(asOnion(x), ...)
}

#' @method summary skymap
#' @export
summary.skymap <- function(object, ...) {
  summary(asOnion(object), ...)
}
