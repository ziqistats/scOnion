## Slots that actually exist on a fitted Onion. The pre-package defaults also
## listed "skymap.feature", a combined gene+peak matrix that Onion() stopped
## returning, so a bare OnionUMAP(object) or OnionFindNeighbors(object) call,
## and the documented slot = "all", always failed on umap(X = NULL).
.onion_slots <- function() {
  c("W", "U", "V", "P", "Q")
}

#' @noRd
.resolve_slots <- function(slot, object) {
  valid <- .onion_slots()
  if (identical(slot, "all")) slot <- valid
  ## POLARIS-era slot names are accepted, with the once-per-session warning.
  leg <- slot %in% names(.legacy_slots)
  if (any(leg)) {
    for (s in slot[leg]) .warn_legacy(s, .legacy_slots[[s]])
    slot[leg] <- unname(.legacy_slots[slot[leg]])
  }
  bad <- setdiff(slot, valid)
  if ("skymap.feature" %in% bad)
    stop(paste0("The 'skymap.feature' slot no longer exists. It was the stacked ",
                "gene-and-feature map, rbind(P * sqrt(nrow(P)), Q * sqrt(nrow(Q))), ",
                "and its construction was commented out of Onion() before this ",
                "code became a package, so no recent fit contains it. Use ",
                "'P' and 'Q' separately, or OnionUMAPChr() / ",
                "OnionFindNeighborsChr() for the per-chromosome maps, which do ",
                "stack the two. If you need the combined map restored as a slot, ",
                "that is a deliberate API decision, not an oversight."),
         call. = FALSE)
  if (length(bad))
    stop(sprintf("Unknown slot(s): %s. Choose from %s, or \"all\".",
                 paste(bad, collapse = ", "), paste(valid, collapse = ", ")),
         call. = FALSE)
  absent <- slot[vapply(slot, function(s) is.null(.subset2(object, s)), logical(1))]
  if (length(absent))
    stop(sprintf("This fit has no %s.", paste(absent, collapse = ", ")),
         call. = FALSE)
  stats::setNames(slot, slot)
}

#' UMAP of an Onion embedding
#'
#' @param object An [Onion][Onion-class] fit from [Onion()].
#' @param slot Which embedding(s) to reduce. One or more of `"W"`,
#'   `"U"`, `"V"`, `"P"`, `"Q"`, or `"all"`.
#'   Defaults to `"W"`, the joint cell embedding.
#' @param n_components,n_neighbors,metric,min_dist Passed to [uwot::umap()].
#'   Note `metric` here takes uwot's lowercase spelling, e.g. `"cosine"`.
#' @param verbose Report progress.
#' @param ... Further arguments for [uwot::umap()].
#'
#' @return The fit with results added under `fit$UMAPs`, one
#'   `data.frame` per slot named `<slot>.UMAP`. Repeated calls accumulate rather
#'   than replace, which is the storage contract the analysis scripts expect
#'   (for example `fit$UMAPs$W.UMAP`).
#'
#' @examples
#' set.seed(1)
#' sim <- onionSimulate(n = 150, g = 50, p = 60, r = 4)
#' fit <- Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE)
#' fit <- OnionUMAP(fit, n_neighbors = 15, verbose = FALSE)
#' head(fit$UMAPs$W.UMAP)
#'
#' @export
OnionUMAP <- function(object,
                       slot = "W",
                       n_components = 3L, n_neighbors = 30L,
                       metric = "cosine", min_dist = 0.2,
                       verbose = TRUE, ...) {

  object <- .check_Onion(object)
  slot <- .resolve_slots(slot, object)

  new <- lapply(slot, function(s) {
    if (isTRUE(verbose)) message("UMAP on ", s, ".")
    m <- as.matrix(object[[s]])
    nn <- min(n_neighbors, nrow(m) - 1L)
    if (nn < 2L)
      stop(sprintf("%s has only %d rows, too few for a UMAP.", s, nrow(m)),
           call. = FALSE)
    co <- as.data.frame(uwot::umap(X = m, n_components = n_components,
                                   n_neighbors = nn, metric = metric,
                                   min_dist = min_dist, ...))
    rownames(co) <- rownames(m)
    colnames(co) <- paste0("UMAP_", seq_len(ncol(co)))
    co
  })
  names(new) <- paste0(slot, ".UMAP")
  object$UMAPs <- .merge_slot(object$UMAPs, new)
  object
}

#' UMAP of the per-chromosome feature embedding
#'
#' Reduces the gene and second-modality loadings of one or more chromosomes,
#' stacking the two feature types into a single UMAP per chromosome so genes and
#' peaks share a coordinate system.
#'
#' @param object An [Onion][Onion-class] fit made with `gene.chr.ref`.
#' @param chr Chromosomes to reduce, or `"all"`.
#' @param n_components,n_neighbors,metric,min_dist Passed to [uwot::umap()].
#' @param verbose Report progress.
#' @param ... Further arguments for [uwot::umap()].
#'
#' @return The fit with `fit$UMAPs$per.chr.UMAP`, a list
#'   keyed by chromosome, each entry a `data.frame` carrying a `feature_type`
#'   column marking rows as `gene` or `feature`.
#'
#' @details
#' This function previously assumed `per.chr[[chr]]` was a single
#' combined matrix. `Onion()` stores `list(P, Q, stdev)` per chromosome, so
#' the loadings are stacked here explicitly.
#'
#' @export
OnionUMAPChr <- function(object, chr = "all",
                           n_components = 2L, n_neighbors = 30L,
                           metric = "cosine", min_dist = 0.2,
                           verbose = TRUE, ...) {

  object <- .check_Onion(object)
  feat <- .resolve_chrs(object, chr)

  out <- lapply(stats::setNames(names(feat), names(feat)), function(ch) {
    if (isTRUE(verbose)) message("UMAP on ", ch, ".")
    cc <- feat[[ch]]
    m <- rbind(as.matrix(cc$P), as.matrix(cc$Q))
    type <- rep(c("gene", "feature"), c(nrow(cc$P), nrow(cc$Q)))
    nn <- min(n_neighbors, nrow(m) - 1L)
    if (nn < 2L) return(NULL)
    co <- as.data.frame(uwot::umap(X = m, n_components = n_components,
                                   n_neighbors = nn, metric = metric,
                                   min_dist = min_dist, ...))
    rownames(co) <- rownames(m)
    colnames(co) <- paste0("UMAP_", seq_len(ncol(co)))
    co$feature_type <- type
    co
  })
  object$UMAPs[["per.chr.UMAP"]] <-
    .merge_slot(object$UMAPs[["per.chr.UMAP"]], out)
  object
}

#' Accumulate named results into an existing slot list
#'
#' Repeated calls to the UMAP and neighbor verbs add to what is already stored
#' rather than replacing it, with a same-named entry overwritten. This is the
#' contract the pre-package code established and that the analysis scripts read.
#' @noRd
.merge_slot <- function(existing, new) {
  if (is.null(existing)) return(new)
  keep <- setdiff(names(existing), names(new))
  c(existing[keep], new)
}

#' @noRd
.resolve_chrs <- function(object, chr) {
  feat <- object$per.chr
  if (is.null(feat))
    stop(paste0("This fit has no per-chromosome feature embedding. Re-run ",
                "Onion() with gene.chr.ref supplied."), call. = FALSE)
  ok <- names(feat)[!vapply(feat, is.null, logical(1))]
  if (!length(ok))
    stop("No chromosome in this fit has a feature embedding.", call. = FALSE)
  if (identical(chr, "all")) chr <- ok
  bad <- setdiff(chr, ok)
  if (length(bad))
    stop(sprintf("No feature embedding for %s. Available: %s.",
                 paste(bad, collapse = ", "), paste(ok, collapse = ", ")),
         call. = FALSE)
  feat[chr]
}
