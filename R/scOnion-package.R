#' scOnion: concordance-aware integration of paired single-cell multiomic data
#'
#' scOnion takes two matrices measured on the same cells and returns, from one
#' closed-form spectral decomposition, a per-cell cross-modality concordance
#' score, a joint cell embedding, and a cross-modal feature embedding that
#' scores candidate feature pairs with an uncertainty.
#'
#' @section Entry point:
#' [Onion()] performs the fit and returns an [Onion][Onion-class] object. Every other
#' function in the package consumes that object:
#'
#' \describe{
#'   \item{[OnionConcordance()]}{per-cell concordance score}
#'   \item{[OnionUMAP()], [OnionFindNeighbors()]}{embeddings and graphs}
#'   \item{[OnionLinkageTable()]}{cis gene-peak linkage with jackknife SE}
#'   \item{[OnionZdev()]}{genomic-distance-aware z-statistic}
#' }
#'
#' @section Notation:
#' The code follows the notation of the scOnion manuscript. \eqn{X} (`n` cells
#' by `g` genes) and \eqn{Y} (`n` cells by `p` features) are the centered and
#' scaled modality matrices. \eqn{\hat{U}(r)} and \eqn{\hat{V}(r)} are the
#' paired cell embeddings, \eqn{W = [\hat{U} | \hat{V}]} the joint cell
#' embedding, and \eqn{P}, \eqn{Q}, \eqn{d} the gene loadings, second-modality
#' loadings and singular-value weights of the joint feature embedding.
#'
#' @section The name:
#' The two denoised modalities are the two halves of an onion, each a nested
#' set of rings, one ring per singular component. The halves meet at a small
#' core whose SVD gives the joint fit, the cross-section \eqn{W} lays the two
#' halves side by side, and the concordance score is the symmetry of the two
#' halves at one cell. The package was called POLARIS until version 0.99.1;
#' the old names still work, see [scOnion-deprecated] and [asOnion()].
#'
#' @keywords internal
#' @useDynLib scOnion, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
