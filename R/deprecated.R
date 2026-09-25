#' Deprecated POLARIS-era names
#'
#' Until version 0.99.1 this package was called POLARIS, its fit was a
#' `skymap`, and its functions were named accordingly. The old names below keep
#' working for one release cycle: each warns that it is deprecated and forwards
#' to its replacement. They will be removed in version 1.0.0. Saved fits in the
#' old layout are handled separately, by [asOnion()], and the old element names
#' of a fit are listed under [Onion-class].
#'
#' | Deprecated | Replacement |
#' |---|---|
#' | `Polaris()` | [Onion()] |
#' | `RunPolaris()` | [RunOnion()]; still returns `list(skymap, object)` and still names the reduction `"polaris"` |
#' | `SkymapSimScore()` | [OnionConcordance()] |
#' | `SkymapLinkageTable()` | [OnionLinkageTable()] |
#' | `SkymapZdev()` | [OnionZdev()] |
#' | `SkymapUMAP()`, `SkymapUMAP.Chr()` | [OnionUMAP()], [OnionUMAPChr()] |
#' | `SkymapFindNeighbors()`, `SkymapFindNeighbors.Chr()` | [OnionFindNeighbors()], [OnionFindNeighborsChr()] |
#' | `polarisTSSRef()`, `polarisMakeTSSRef()` | [onionTSSRef()], [onionMakeTSSRef()] |
#' | `polarisSimulate()` | [onionSimulate()] |
#' | `polarisFitSummary()`, `polarisVersions()`, `polarisSpectrum()` | [onionFitSummary()], [onionVersions()], [onionSpectrum()] |
#' | `POLARIS_DIST_BREAKS`, `POLARIS_MIN_DIST` | [ONION_DIST_BREAKS], [ONION_MIN_DIST] |
#'
#' @param ... Arguments passed to the replacement.
#' @param skymap,tab,object The first argument of the replacement.
#' @param reduction.name Kept at `"polaris"` here, so that code written against
#'   `RunPolaris()` keeps finding its reduction. [RunOnion()] uses `"onion"`.
#'
#' @name scOnion-deprecated
#' @keywords internal
#' @include zdev.R
NULL

#' @rdname scOnion-deprecated
#' @export
Polaris <- function(...) {
  .Deprecated("Onion", package = "scOnion")
  Onion(...)
}

#' @rdname scOnion-deprecated
#' @export
RunPolaris <- function(object, ..., reduction.name = "polaris") {
  .Deprecated("RunOnion", package = "scOnion")
  r <- RunOnion(object, ..., reduction.name = reduction.name)
  list(skymap = r$fit, object = r$object)
}

#' @rdname scOnion-deprecated
#' @export
SkymapSimScore <- function(skymap, ...) {
  .Deprecated("OnionConcordance", package = "scOnion")
  OnionConcordance(skymap, ...)
}

#' @rdname scOnion-deprecated
#' @export
SkymapLinkageTable <- function(skymap, ...) {
  .Deprecated("OnionLinkageTable", package = "scOnion")
  OnionLinkageTable(skymap, ...)
}

#' @rdname scOnion-deprecated
#' @export
SkymapZdev <- function(tab, ...) {
  .Deprecated("OnionZdev", package = "scOnion")
  OnionZdev(tab, ...)
}

#' @rdname scOnion-deprecated
#' @export
SkymapUMAP <- function(skymap, ...) {
  .Deprecated("OnionUMAP", package = "scOnion")
  OnionUMAP(skymap, ...)
}

#' @rdname scOnion-deprecated
#' @export
SkymapUMAP.Chr <- function(skymap, ...) {
  .Deprecated("OnionUMAPChr", package = "scOnion")
  OnionUMAPChr(skymap, ...)
}

#' @rdname scOnion-deprecated
#' @export
SkymapFindNeighbors <- function(skymap, ...) {
  .Deprecated("OnionFindNeighbors", package = "scOnion")
  OnionFindNeighbors(skymap, ...)
}

#' @rdname scOnion-deprecated
#' @export
SkymapFindNeighbors.Chr <- function(skymap, ...) {
  .Deprecated("OnionFindNeighborsChr", package = "scOnion")
  OnionFindNeighborsChr(skymap, ...)
}

#' @rdname scOnion-deprecated
#' @export
polarisTSSRef <- function(...) {
  .Deprecated("onionTSSRef", package = "scOnion")
  onionTSSRef(...)
}

#' @rdname scOnion-deprecated
#' @export
polarisMakeTSSRef <- function(...) {
  .Deprecated("onionMakeTSSRef", package = "scOnion")
  onionMakeTSSRef(...)
}

#' @rdname scOnion-deprecated
#' @export
polarisSimulate <- function(...) {
  .Deprecated("onionSimulate", package = "scOnion")
  onionSimulate(...)
}

#' @rdname scOnion-deprecated
#' @export
polarisFitSummary <- function(...) {
  .Deprecated("onionFitSummary", package = "scOnion")
  onionFitSummary(...)
}

#' @rdname scOnion-deprecated
#' @export
polarisVersions <- function(...) {
  .Deprecated("onionVersions", package = "scOnion")
  onionVersions(...)
}

#' @rdname scOnion-deprecated
#' @export
polarisSpectrum <- function(skymap, ...) {
  .Deprecated("onionSpectrum", package = "scOnion")
  onionSpectrum(skymap, ...)
}

#' @rdname scOnion-deprecated
#' @format NULL
#' @export
POLARIS_DIST_BREAKS <- ONION_DIST_BREAKS

#' @rdname scOnion-deprecated
#' @format NULL
#' @export
POLARIS_MIN_DIST <- ONION_MIN_DIST
