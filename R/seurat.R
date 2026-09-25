#' Fit scOnion to a Seurat object
#'
#' Convenience wrapper that pulls the two modality matrices out of a Seurat
#' object, orients them as scOnion expects, and calls [Onion()]. It also
#' returns the joint cell embedding as a Seurat `DimReduc`, so the fit can be
#' used directly with `FindNeighbors`, `FindClusters` and `RunUMAP`.
#'
#' `Onion()` itself takes plain matrices and is the reproducible entry point.
#' This function only assembles its arguments.
#'
#' @param object A `Seurat` object with both modalities as assays.
#' @param assay.x,assay.y Assay names. Defaults `"RNA"` and `"ATAC"`.
#' @param layer.x,layer.y Layer (slot) to take from each assay. The defaults are
#'   the conventions used throughout the scOnion manuscript: `"scale.data"` for
#'   RNA, so \eqn{X} is centered and scaled, and `"data"` for ATAC, which is
#'   expected to hold TF-IDF normalized values from `Signac::RunTFIDF`. For
#'   CITE-seq use `assay.y = "ADT"` with CLR-normalized `"data"`.
#' @param gene.chr.ref Gene annotation. If `NULL` and `assay.y` carries a Signac
#'   annotation, that is used, which is the safest source because it is
#'   guaranteed to match the genome the peaks were called on. Pass a `GRanges`
#'   to override, or `FALSE` to skip the per-chromosome step.
#' @param reduction.name Name under which to store the joint cell embedding.
#'   Its key is `onion_`, matching the column names of `fit$W`.
#' @param ... Passed to [Onion()], including `x.sds`, `T1`, `T2`, `T3`,
#'   `mc.cores` and `seed`.
#'
#' @return A list with `fit` (the [Onion][Onion-class] fit) and `object` (the input
#'   Seurat object with the joint embedding added as `reduction.name`).
#'
#' @section Preprocessing:
#' scOnion does not normalize for you, because the normalization is a
#' methodological choice that belongs to your analysis. The conventions used in
#' the manuscript are:
#'
#' ```
#' obj <- NormalizeData(obj, assay = "RNA",
#'                      normalization.method = "LogNormalize",
#'                      scale.factor = 1e6)
#' obj <- FindVariableFeatures(obj, nfeatures = 3000)
#' # NOTE: ScaleData is deliberately run on ALL genes, not on the variable
#' # features. scOnion selects its own feature space by the x.sds variance
#' # floor, which is broader than a highly-variable-gene list.
#' obj <- ScaleData(obj, assay = "RNA", features = rownames(obj))
#' DefaultAssay(obj) <- "ATAC"
#' obj <- Signac::RunTFIDF(obj)
#' ```
#'
#' Note the interaction with `x.sds`: its default of `0.95` assumes \eqn{X} is
#' z-scored. If you pass a log-normalized layer instead, nearly every gene falls
#' below the floor and the fit will stop with an empty \eqn{X}. Pass
#' `x.sds = 0` in that case.
#'
#' @examples
#' \dontrun{
#' res <- RunOnion(obj, assay.x = "RNA", assay.y = "ATAC", mc.cores = 4)
#' res$object <- Seurat::FindNeighbors(res$object, reduction = "onion",
#'                                     dims = 1:ncol(res$fit$W))
#' res$fit <- OnionConcordance(res$fit)
#' res$object$concordance <- res$fit$concordance[colnames(res$object)]
#' }
#'
#' @seealso [Onion()]
#' @export
RunOnion <- function(object,
                       assay.x = "RNA", assay.y = "ATAC",
                       layer.x = "scale.data", layer.y = "data",
                       gene.chr.ref = NULL,
                       reduction.name = "onion",
                       ...) {

  for (p in c("Seurat", "SeuratObject"))
    if (!requireNamespace(p, quietly = TRUE))
      stop(sprintf("RunOnion() needs the %s package.", p), call. = FALSE)
  if (!methods::is(object, "Seurat"))
    stop("`object` must be a Seurat object. For matrices use Onion().",
         call. = FALSE)

  for (a in c(assay.x, assay.y))
    if (!a %in% names(object@assays))
      stop(sprintf("No assay named '%s'. Available: %s.", a,
                   paste(names(object@assays), collapse = ", ")), call. = FALSE)

  grab <- function(assay, layer) {
    m <- SeuratObject::LayerData(object, assay = assay, layer = layer)
    if (is.null(m) || !length(m))
      stop(sprintf(paste0("Layer '%s' of assay '%s' is empty. For scale.data, ",
                          "run ScaleData(obj, assay = '%s', features = ",
                          "rownames(obj)) first."), layer, assay, assay),
           call. = FALSE)
    Matrix::t(m)                      # features x cells -> cells x features
  }

  X <- grab(assay.x, layer.x)
  Y <- grab(assay.y, layer.y)
  ## Seurat can hold different cell sets per assay; align on the intersection in
  ## the object's own column order.
  common <- intersect(rownames(X), rownames(Y))
  if (!length(common))
    stop("The two assays share no cells.", call. = FALSE)
  if (length(common) < nrow(X) || length(common) < nrow(Y))
    message(sprintf("Using the %d cells present in both assays.", length(common)))
  X <- X[common, , drop = FALSE]
  Y <- Y[common, , drop = FALSE]

  if (isFALSE(gene.chr.ref)) {
    gene.chr.ref <- NULL
  } else if (is.null(gene.chr.ref)) {
    gene.chr.ref <- tryCatch({
      a <- object[[assay.y]]
      if (!is.null(methods::slot(a, "annotation"))) methods::slot(a, "annotation") else NULL
    }, error = function(e) NULL)
    if (!is.null(gene.chr.ref)) {
      message("Using the annotation attached to assay '", assay.y, "'.")
      if (!"gene_name" %in% names(S4Vectors::mcols(gene.chr.ref)))
        stop(paste0("The attached annotation has no `gene_name` column. Pass ",
                    "gene.chr.ref explicitly, e.g. onionTSSRef('hg38')."),
             call. = FALSE)
    } else {
      message(paste0("No annotation attached to assay '", assay.y, "', so the ",
                     "per-chromosome step is skipped. Pass gene.chr.ref to ",
                     "enable OnionLinkageTable()."))
    }
  }

  fit <- Onion(X, Y, gene.chr.ref = gene.chr.ref, ...)

  object[[reduction.name]] <- SeuratObject::CreateDimReducObject(
    embeddings = fit$W,
    stdev      = as.numeric(fit$stdev.cell),
    key        = "onion_",
    assay      = assay.x)

  list(fit = fit, object = object)
}
