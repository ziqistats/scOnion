## Each test below that references a bug is pinning behaviour that was wrong in
## the pre-package scripts, so a regression would be caught rather than
## rediscovered.

sim_fit <- function(...) {
  set.seed(11)
  sim <- onionSimulate(n = 150, g = 50, p = 70, r = 4)
  list(sim = sim, fit = Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE, ...))
}

test_that("Onion returns a well-formed object", {
  z <- sim_fit()
  fit <- z$fit
  expect_s3_class(fit, "Onion")
  expect_equal(nrow(fit$W), 150)
  expect_equal(ncol(fit$W), ncol(fit$U) + ncol(fit$V))
  expect_identical(rownames(fit$W), rownames(z$sim$X))
  expect_identical(rownames(fit$P), colnames(z$sim$X))
  expect_identical(rownames(fit$Q), colnames(z$sim$Y))
  ## U, V, P, Q must be matrices, not data.frames: as data.frames,
  ## OnionConcordance(metric = "cor") silently produced an all-NA result of the
  ## wrong length.
  for (s in c("U", "V", "P", "Q"))
    expect_true(is.matrix(fit[[s]]), info = s)
})

test_that("the fit records the ranks it actually used", {
  fit <- sim_fit()$fit
  ip <- fit$input.param
  ## The pre-package version captured input.param before any default resolved,
  ## so the saved object recorded only the thresholds and the ranks were NULL.
  for (nm in c("r1", "r2", "r3", "r4", "r5"))
    expect_true(is.numeric(ip[[nm]]) && length(ip[[nm]]) == 1L, info = nm)
  expect_equal(ip$r3, ncol(fit$U))
  expect_equal(ip$r5, ncol(fit$P))
})

test_that("Onion is deterministic and leaves the caller's RNG alone", {
  set.seed(11); sim <- onionSimulate(n = 120, g = 40, p = 50, r = 3)
  a <- Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE)
  b <- Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE)
  expect_equal(a$stdev.cell, b$stdev.cell, tolerance = 0)

  ## irlba draws a random start, so without the internal seed handling a fit
  ## would advance the user's stream and break their reproducibility.
  set.seed(99); before <- stats::runif(1)
  set.seed(99); invisible(Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE))
  after <- stats::runif(1)
  expect_equal(before, after)
})

test_that("mismatched or unnamed cells are caught", {
  set.seed(11); sim <- onionSimulate(n = 80, g = 30, p = 40, r = 3)
  expect_error(Onion(sim$X, sim$Y[1:79, ], verbose = FALSE),
               "same number of cells")
  shuffled <- sim$Y[sample(nrow(sim$Y)), ]
  expect_error(Onion(sim$X, shuffled, x.sds = 0, verbose = FALSE),
               "differently ordered cell names")
  bare.x <- sim$X; bare.y <- sim$Y
  rownames(bare.x) <- rownames(bare.y) <- NULL
  expect_warning(Onion(bare.x, bare.y, x.sds = 0, verbose = FALSE),
                 "no cell names")
})

test_that("an over-strict variance floor is reported, not left to fail later", {
  set.seed(11); sim <- onionSimulate(n = 80, g = 30, p = 40, r = 3)
  ## x.sds = 0.95 on non-z-scored input silently emptied X in the pre-package
  ## code, and the failure surfaced much further downstream.
  expect_error(Onion(sim$X, sim$Y, x.sds = 100, verbose = FALSE),
               "standard deviation above x.sds")
})

test_that("a flat spectrum gives a named error, not '1:-Inf'", {
  ## max(which(...)) returned -Inf when no ratio cleared the threshold, and the
  ## fit then died on `1:r` with "result would be too long a vector".
  set.seed(11); sim <- onionSimulate(n = 100, g = 40, p = 50, r = 3)
  expect_error(Onion(sim$X, sim$Y, T1 = 1e6, x.sds = 0, verbose = FALSE),
               "No spectral gap")
})

test_that("modalities with 50 or fewer features still fit", {
  ## min(50, ncol(X)) equalled min(dim) here, which made irlba hard-error.
  set.seed(11); sim <- onionSimulate(n = 200, g = 40, p = 45, r = 4)
  expect_no_error(Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE))
})

test_that("gene.chr.ref is optional and its absence is explained", {
  set.seed(11); sim <- onionSimulate(n = 120, g = 40, p = 50, r = 3)
  fit <- Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE)
  expect_null(fit$per.chr)
  ## The default gene.chr.ref = NULL used to fail inside mclapply and return a
  ## list of try-error objects that looked like a successful fit.
  expect_error(OnionLinkageTable(fit, gene.anno = sim$anno, verbose = FALSE),
               "no per-chromosome feature embedding")
})

test_that("the per-chromosome embedding is built when an annotation is given", {
  set.seed(11); sim <- onionSimulate(n = 200, g = 60, p = 90, r = 5)
  fit <- Onion(sim$X, sim$Y, gene.chr.ref = sim$anno, x.sds = 0, verbose = FALSE)
  expect_true(is.list(fit$per.chr))
  expect_true("chr1" %in% names(fit$per.chr))
  cc <- fit$per.chr[["chr1"]]
  expect_named(cc, c("P", "Q", "stdev"))
  expect_equal(nrow(cc$P) , 60)
})

test_that("the retired skymap.feature slot fails with an explanation", {
  ## Six analysis scripts still reference this slot, so the error has to say
  ## what happened rather than just "unknown slot".
  fit <- sim_fit()$fit
  expect_error(OnionUMAP(fit, slot = "skymap.feature"), "no longer exists")
  expect_error(OnionFindNeighbors(fit, slot = "skymap.feature"), "no longer exists")
  expect_error(OnionUMAP(fit, slot = "nonsense"), "Unknown slot")
})

test_that("r3 and r4 are validated against the available rank", {
  set.seed(11); sim <- onionSimulate(n = 150, g = 50, p = 60, r = 4)
  expect_error(Onion(sim$X, sim$Y, r3 = 500, x.sds = 0, verbose = FALSE),
               "cannot exceed")
  expect_error(Onion(sim$X, sim$Y, r3 = 0, x.sds = 0, verbose = FALSE),
               "at least 1")
})

test_that("n.pc is clamped against both U and V when r3 != r4", {
  set.seed(11); sim <- onionSimulate(n = 200, g = 60, p = 70, r = 5)
  fit <- Onion(sim$X, sim$Y, r3 = 5, r4 = 3, x.sds = 0, verbose = FALSE)
  expect_equal(ncol(fit$U), 5)
  expect_equal(ncol(fit$V), 3)
  ## slicing V by ncol(U) would be out of bounds
  out <- OnionConcordance(fit)
  expect_length(out$concordance, 200)
  expect_false(anyNA(out$concordance))
})

test_that("UMAP and neighbor results use the $UMAPs / $NN.graph contract", {
  ## The analysis scripts read object$UMAPs$<slot>.UMAP and
  ## object$NN.graph$<slot>.nb (PLOT_motif_final.R:136, PLOT_HVG.R:43,
  ## PLOT4_nbh_distance.R:46, validation_linkage_methods.R:137). A flat naming
  ## scheme would silently break all of them.
  set.seed(11); sim <- onionSimulate(n = 150, g = 50, p = 60, r = 4)
  fit <- Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE)

  fit <- OnionUMAP(fit, n_neighbors = 15, verbose = FALSE)
  expect_true(is.list(fit$UMAPs))
  expect_true("W.UMAP" %in% names(fit$UMAPs))
  expect_equal(nrow(fit$UMAPs$W.UMAP), 150)

  fit <- OnionFindNeighbors(fit, n.neighbors = 10, verbose = FALSE)
  expect_true(is.list(fit$NN.graph))
  expect_true("W.nb" %in% names(fit$NN.graph))
  expect_equal(nrow(fit$NN.graph$W.nb$index), 150)

  ## repeated calls must ACCUMULATE, not replace
  fit <- OnionUMAP(fit, slot = "U", n_neighbors = 15, verbose = FALSE)
  expect_true(all(c("W.UMAP", "U.UMAP") %in% names(fit$UMAPs)))
  fit <- OnionFindNeighbors(fit, slot = "U", n.neighbors = 10,
                             verbose = FALSE)
  expect_true(all(c("W.nb", "U.nb") %in% names(fit$NN.graph)))
})

test_that("the per-chromosome verbs use the nested contract too", {
  set.seed(11); sim <- onionSimulate(n = 200, g = 60, p = 90, r = 5)
  fit <- Onion(sim$X, sim$Y, gene.chr.ref = sim$anno, x.sds = 0,
                 verbose = FALSE)
  fit <- OnionFindNeighborsChr(fit, n.neighbors = 10, verbose = FALSE)
  expect_true("per.chr.nb" %in% names(fit$NN.graph))
  expect_true("chr1" %in% names(fit$NN.graph$per.chr.nb))
  expect_named(fit$NN.graph$per.chr.nb$chr1,
               c("gene.to.feature", "feature.to.gene"))
})

test_that("reconcile.ranks controls the rank collapse", {
  ## The scOnion manuscript mixes conventions: Figs 2, 3, 4b,c and 5 were fitted
  ## WITHOUT the collapse, Fig 6 and the per-cell-type fits WITH it. The package
  ## must be able to reproduce both, so this is an argument rather than a
  ## hard-coded choice (measured audit, 2026-09-10/11).
  set.seed(11); sim <- onionSimulate(n = 250, g = 70, p = 90, r = 6)

  a <- Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE, reconcile.ranks = TRUE)
  b <- Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE, reconcile.ranks = FALSE)

  expect_true(a$input.param$reconcile.ranks)
  expect_false(b$input.param$reconcile.ranks)
  ## collapsed: the two modality ranks are equal. unreconciled: each keeps its own.
  expect_equal(a$input.param$r1, a$input.param$r2)
  expect_equal(b$input.param$r1, b$input.param$r1.auto)
  expect_equal(b$input.param$r2, b$input.param$r2.auto)

  ## r3 = r4 = min(r1, r2) either way, so the joint embedding has the SAME width
  ## under both conventions -- ncol(W) cannot reveal which was used.
  expect_equal(ncol(a$W), ncol(b$W))
})

test_that("an unequal-rank fit actually differs from the collapsed one", {
  ## Only meaningful when the ratio rule picks different ranks for the two
  ## modalities; construct a case where it does.
  set.seed(23); sim <- onionSimulate(n = 300, g = 80, p = 100, r = 6)
  b <- Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE, reconcile.ranks = FALSE)
  skip_if(b$input.param$r1 == b$input.param$r2,
          "ratio rule happened to pick equal ranks; nothing to contrast")
  a <- Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE, reconcile.ranks = TRUE)
  ## the embeddings must genuinely differ, else the argument is a no-op
  expect_false(isTRUE(all.equal(as.numeric(a$stdev.cell),
                                as.numeric(b$stdev.cell))))
})

test_that("onionFitSummary records which convention was used", {
  ## Supplementary Table 3 has to state this per analysis, and it is invisible
  ## from the embedding width.
  set.seed(11); sim <- onionSimulate(n = 200, g = 60, p = 80, r = 5)
  a <- Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE, reconcile.ranks = TRUE)
  b <- Onion(sim$X, sim$Y, x.sds = 0, verbose = FALSE, reconcile.ranks = FALSE)
  s <- onionFitSummary(collapsed = a, unreconciled = b)
  expect_true("reconcile.ranks" %in% names(s))
  expect_equal(s$reconcile.ranks, c(TRUE, FALSE))
})
