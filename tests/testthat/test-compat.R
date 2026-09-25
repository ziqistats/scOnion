## The package was POLARIS until 0.99.1. Terabytes of fits are saved under the
## old element names and hundreds of scripts call the old functions, so the
## rename ships with a compatibility layer for one release cycle. These tests
## pin it: saved fits upgrade, old names still read, old functions forward.

compat_fit <- function() {
  set.seed(5)
  sim <- onionSimulate(n = 150, g = 50, p = 60, r = 4)
  list(sim = sim,
       fit = Onion(sim$X, sim$Y, gene.chr.ref = sim$anno, x.sds = 0,
                   verbose = FALSE))
}

## Exactly the layout POLARIS 0.99.1 saved: class "skymap", old element names,
## and the version under polaris.version.
as_legacy <- function(fit) {
  old <- unclass(fit)
  map <- c(W = "skymap.cell", U = "skymap.U", V = "skymap.V",
           P = "skymap.P", Q = "skymap.Q", per.chr = "skymap.feature.chr")
  nm <- names(old)
  hit <- nm %in% names(map)
  nm[hit] <- unname(map[nm[hit]])
  names(old) <- nm
  old$input.param$polaris.version <- old$input.param$onion.version
  old$input.param$onion.version <- NULL
  class(old) <- c("skymap", "list")
  old
}

## The legacy warnings and message fire once per session; reset between tests.
reset_legacy <- function() {
  e <- scOnion:::.legacy_env
  rm(list = ls(envir = e), envir = e)
}

test_that("asOnion() restores the current layout from a 0.99.1 fit", {
  z <- compat_fit()
  old <- as_legacy(z$fit)
  expect_s3_class(old, "skymap")
  expect_null(old$W)

  new <- asOnion(old)
  expect_s3_class(new, "Onion")
  for (s in c("W", "U", "V", "P", "Q", "per.chr", "stdev.cell",
              "stdev.feature", "X.svd", "Y.svd"))
    expect_identical(new[[s]], z$fit[[s]], info = s)
  expect_identical(new$input.param$onion.version,
                   z$fit$input.param$onion.version)
  ## the record of what produced the fit is kept
  expect_identical(new$input.param$polaris.version,
                   z$fit$input.param$onion.version)
  ## idempotent, and the pre-package bare list converts too
  expect_identical(asOnion(new), new)
  expect_s3_class(asOnion(unclass(old)), "Onion")
  expect_error(asOnion(list(a = 1)), "must be an Onion fit")
  expect_error(asOnion(42), "must be an Onion fit")
})

test_that("asOnion() renames the derived UMAP and neighbor keys", {
  z <- compat_fit()
  fit <- OnionUMAP(z$fit, n_neighbors = 15, verbose = FALSE)
  fit <- OnionFindNeighbors(fit, n.neighbors = 10, verbose = FALSE)
  fit <- OnionFindNeighborsChr(fit, n.neighbors = 5, verbose = FALSE)
  old <- as_legacy(fit)
  names(old$UMAPs) <- "skymap.cell.UMAP"
  names(old$NN.graph) <- c("skymap.cell.nb", "skymap.feature.chr.nb")
  new <- asOnion(old)
  expect_named(new$UMAPs, "W.UMAP")
  expect_named(new$NN.graph, c("W.nb", "per.chr.nb"))
})

test_that("every verb accepts a legacy fit and upgrades it on entry", {
  z <- compat_fit()
  old <- as_legacy(z$fit)

  reset_legacy()
  expect_message(sc <- OnionConcordance(old), "asOnion")
  expect_s3_class(sc, "Onion")
  expect_equal(sc$concordance, OnionConcordance(z$fit)$concordance)
  ## the message is once per session
  expect_silent(OnionConcordance(old))

  tab.old <- OnionLinkageTable(old, gene.anno = z$sim$anno, verbose = FALSE)
  tab.new <- OnionLinkageTable(z$fit, gene.anno = z$sim$anno, verbose = FALSE)
  expect_equal(tab.old, tab.new)

  expect_equal(onionSpectrum(old), onionSpectrum(z$fit))
  expect_true("W.UMAP" %in%
                names(OnionUMAP(old, n_neighbors = 15, verbose = FALSE)$UMAPs))
  expect_true("W.nb" %in%
                names(OnionFindNeighbors(old, n.neighbors = 10,
                                         verbose = FALSE)$NN.graph))
  expect_output(print(old), "An Onion")
  expect_output(summary(old), "An Onion")
})

test_that("onionFitSummary reads the version from a POLARIS-era fit", {
  z <- compat_fit()
  old <- as_legacy(z$fit)
  s <- onionFitSummary(old = old)
  expect_true("onion.version" %in% names(s))
  expect_identical(s$onion.version, z$fit$input.param$onion.version)
  ## a single legacy fit is one fit, not a list of fits
  expect_equal(nrow(onionFitSummary(old)), 1)
})

test_that("legacy element names still read through $ and [[, warning once", {
  z <- compat_fit()
  fit <- z$fit

  reset_legacy()
  expect_warning(w <- fit$skymap.cell, "POLARIS-era name")
  expect_identical(w, fit$W)
  expect_silent(w2 <- fit$skymap.cell)          # second access is quiet
  expect_identical(w2, fit$W)

  reset_legacy()
  expect_warning(p <- fit[["skymap.P"]], "POLARIS-era name")
  expect_identical(p, fit$P)
  expect_warning(pc <- fit$skymap.feature.chr, "POLARIS-era name")
  expect_identical(pc, fit$per.chr)

  ## ordinary access is untouched
  expect_null(fit$nonsense)
  expect_identical(fit[["stdev.cell"]], fit$stdev.cell)
  expect_identical(fit[[1L]], fit$W)
})

test_that("legacy slot strings are accepted by slot =", {
  z <- compat_fit()
  reset_legacy()
  expect_warning(fit <- OnionUMAP(z$fit, slot = "skymap.cell",
                                  n_neighbors = 15, verbose = FALSE),
                 "POLARIS-era name")
  expect_true("W.UMAP" %in% names(fit$UMAPs))
  expect_warning(fit <- OnionFindNeighbors(fit, slot = "skymap.U",
                                           n.neighbors = 10, verbose = FALSE),
                 "POLARIS-era name")
  expect_true("U.nb" %in% names(fit$NN.graph))
  ## the slot retired before the rename still gets its explanation
  expect_error(OnionUMAP(z$fit, slot = "skymap.feature"), "no longer exists")
})

test_that("the deprecated POLARIS names forward to their replacements", {
  z <- compat_fit()
  sim <- z$sim

  expect_warning(a <- Polaris(sim$X, sim$Y, gene.chr.ref = sim$anno,
                              x.sds = 0, verbose = FALSE), "deprecated")
  expect_identical(a, z$fit)

  ## the old first-argument name must still be accepted by name
  expect_warning(b <- SkymapSimScore(skymap = z$fit), "deprecated")
  expect_equal(b, OnionConcordance(z$fit))

  expect_warning(t1 <- SkymapLinkageTable(z$fit, gene.anno = sim$anno,
                                          se = FALSE, verbose = FALSE),
                 "deprecated")
  expect_equal(t1, OnionLinkageTable(z$fit, gene.anno = sim$anno, se = FALSE,
                                     verbose = FALSE))

  tab <- data.frame(gene = "A", peak = "chr1-1-2", chr = "chr1",
                    dist = c(2e4, 3e4), score_raw = c(1, 2), SE = c(1, 1))
  expect_warning(zt <- SkymapZdev(tab), "deprecated")
  expect_equal(zt, OnionZdev(tab))

  expect_warning(u <- SkymapUMAP(z$fit, n_neighbors = 15, verbose = FALSE),
                 "deprecated")
  expect_true("W.UMAP" %in% names(u$UMAPs))
  expect_warning(uc <- SkymapUMAP.Chr(z$fit, n_neighbors = 5, verbose = FALSE),
                 "deprecated")
  expect_true("per.chr.UMAP" %in% names(uc$UMAPs))
  expect_warning(nb <- SkymapFindNeighbors(z$fit, n.neighbors = 10,
                                           verbose = FALSE), "deprecated")
  expect_true("W.nb" %in% names(nb$NN.graph))
  expect_warning(nbc <- SkymapFindNeighbors.Chr(z$fit, n.neighbors = 5,
                                                verbose = FALSE), "deprecated")
  expect_true("per.chr.nb" %in% names(nbc$NN.graph))

  expect_warning(s <- polarisSimulate(n = 20, g = 10, p = 12, r = 2),
                 "deprecated")
  expect_named(s, c("X", "Y", "Z", "shared", "anno"))
  expect_warning(fs <- polarisFitSummary(one = z$fit), "deprecated")
  expect_equal(fs, onionFitSummary(one = z$fit))
  expect_warning(v <- polarisVersions(), "deprecated")
  expect_equal(v, onionVersions())
  expect_warning(sp <- polarisSpectrum(z$fit), "deprecated")
  expect_equal(sp, onionSpectrum(z$fit))
  expect_warning(an <- polarisTSSRef("hg38"), "deprecated")
  expect_equal(length(an), length(onionTSSRef("hg38")))

  expect_identical(POLARIS_DIST_BREAKS, ONION_DIST_BREAKS)
  expect_identical(POLARIS_MIN_DIST, ONION_MIN_DIST)
})
