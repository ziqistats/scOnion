#' Nearest-neighbor graphs on Onion embeddings
#'
#' @param object An [Onion][Onion-class] fit from [Onion()].
#' @param slot Which embedding(s) to build a graph on. One or more of
#'   `"W"`, `"U"`, `"V"`, `"P"`, `"Q"`, or
#'   `"all"`. Defaults to `"W"`.
#' @param n.neighbors Neighbors per point.
#' @param metric Distance for the Annoy backend. Note the capitalized spelling
#'   BiocNeighbors expects, e.g. `"Cosine"`, which differs from uwot's
#'   `"cosine"` in [OnionUMAP()].
#' @param mc.cores Cores, one slot per core.
#' @param verbose Report progress.
#' @param ... Further arguments for [BiocNeighbors::findKNN()].
#'
#' @return The fit with results added under `fit$NN.graph`, one entry
#'   per slot named `<slot>.nb`, each carrying `index` and `distance` matrices
#'   with rownames taken from the embedding. Repeated calls accumulate, with a
#'   same-named entry replaced, matching the storage contract the analysis
#'   scripts read (for example `fit$NN.graph$W.nb`).
#'
#' @examples
#' set.seed(1)
#' sim <- onionSimulate(n = 150, g = 50, p = 60, r = 4)
#' fit <- Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE)
#' fit <- OnionFindNeighbors(fit, n.neighbors = 10, verbose = FALSE)
#' dim(fit$NN.graph$W.nb$index)
#'
#' @export
OnionFindNeighbors <- function(object,
                                slot = "W",
                                n.neighbors = 300L,
                                metric = "Cosine",
                                mc.cores = 1L,
                                verbose = TRUE, ...) {

  object <- .check_Onion(object)
  slot <- .resolve_slots(slot, object)

  res <- .maybe_mclapply(slot, function(s) {
    if (isTRUE(verbose)) message("Neighbors on ", s, ".")
    .knn(as.matrix(object[[s]]), n.neighbors, metric, ...)
  }, mc.cores, "neighbor graphs")

  names(res) <- paste0(slot, ".nb")
  object$NN.graph <- .merge_slot(object$NN.graph, res)
  object
}

#' Cross-modality nearest neighbors within each chromosome
#'
#' For every chromosome, finds each gene's nearest second-modality features and
#' each feature's nearest genes in the shared per-chromosome loading space.
#'
#' @param object An [Onion][Onion-class] fit made with `gene.chr.ref`.
#' @param chr Chromosomes to process, or `"all"`.
#' @param n.neighbors Neighbors per query.
#' @param metric Distance for the Annoy backend, e.g. `"Cosine"`.
#' @param mc.cores Cores, one chromosome per core.
#' @param verbose Report progress.
#'
#' @return The fit with `fit$NN.graph$per.chr.nb`, a list
#'   keyed by chromosome, each holding `gene.to.feature` and `feature.to.gene`.
#'
#' @details
#' Rewritten for the `list(P, Q, stdev)` layout that `Onion()` stores. The
#' pre-package version indexed `per.chr[[chr]]` as a single matrix
#' and called `rownames()` on it, which is `NULL` for a list, so it could not
#' run on any current fit.
#'
#' @export
OnionFindNeighborsChr <- function(object, chr = "all",
                                    n.neighbors = 200L,
                                    metric = "Cosine",
                                    mc.cores = 1L,
                                    verbose = TRUE) {

  object <- .check_Onion(object)
  feat <- .resolve_chrs(object, chr)

  res <- .maybe_mclapply(names(feat), function(ch) {
    if (isTRUE(verbose)) message("Cross-modality neighbors on ", ch, ".")
    cc <- feat[[ch]]
    P <- as.matrix(cc$P); Q <- as.matrix(cc$Q)
    kg <- min(n.neighbors, nrow(P))
    kq <- min(n.neighbors, nrow(Q))
    list(
      gene.to.feature = .knn_query(Q, P, kq, metric, rownames(P), rownames(Q)),
      feature.to.gene = .knn_query(P, Q, kg, metric, rownames(Q), rownames(P)))
  }, mc.cores, "chromosomes")

  names(res) <- names(feat)
  object$NN.graph[["per.chr.nb"]] <-
    .merge_slot(object$NN.graph[["per.chr.nb"]], res)
  object
}

#' @noRd
.knn <- function(m, k, metric, ...) {
  k <- min(k, nrow(m) - 1L)
  if (k < 1L)
    stop(sprintf("Need at least 2 rows to build a neighbor graph; got %d.",
                 nrow(m)), call. = FALSE)
  r <- BiocNeighbors::findKNN(m, k = k,
                              BNPARAM = BiocNeighbors::AnnoyParam(distance = metric),
                              ...)
  rownames(r$index) <- rownames(m)
  rownames(r$distance) <- rownames(m)
  r
}

#' @noRd
.knn_query <- function(reference, query, k, metric, query.names, ref.names) {
  k <- min(k, nrow(reference))
  if (k < 1L) return(NULL)
  r <- BiocNeighbors::queryKNN(
    X = reference, query = query, k = k,
    BNPARAM = BiocNeighbors::AnnoyParam(distance = metric))
  rownames(r$index) <- query.names
  rownames(r$distance) <- query.names
  ## Translate reference row numbers into feature names, which is what a caller
  ## actually wants from a cross-modality graph.
  r$name <- matrix(ref.names[r$index], nrow = nrow(r$index),
                   dimnames = list(query.names, NULL))
  r
}
