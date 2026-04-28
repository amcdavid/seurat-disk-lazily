# Tests for LoadH5SCE.
#
# These exercise the new SingleCellExperiment loader: parameter-validation
# rules (only counts + metadata supported), the cells= subset path, and
# (when HDF5Array is installed) the as.delayed = TRUE disk-backed mode.

skip_if_no_sce <- function() {
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")
}

# ---------------------------------------------------------------------------
# Parameter validation
# ---------------------------------------------------------------------------

test_that("LoadH5SCE rejects unsupported components when not FALSE", {
  skip_if_no_sce()
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  for (arg in c("reductions", "graphs", "neighbors", "images",
                "commands", "misc", "tools")) {
    args <- list(file = path, verbose = FALSE)
    args[[arg]] <- TRUE
    expect_error(do.call(LoadH5SCE, args), regexp = arg)
  }
})

test_that("LoadH5SCE rejects non-character file argument", {
  skip_if_no_sce()
  expect_error(LoadH5SCE(123, verbose = FALSE), regexp = "path")
})

# ---------------------------------------------------------------------------
# Counts + metadata roundtrip
# ---------------------------------------------------------------------------

test_that("LoadH5SCE returns full counts and metadata when no subset", {
  skip_if_no_sce()
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  sce <- LoadH5SCE(path, verbose = FALSE)

  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(ncol(sce), N_CELLS)
  expect_equal(nrow(sce), N_FEATURES)
  expect_equal(colnames(sce), colnames(make_test_seurat()))
  expect_true("counts" %in% SummarizedExperiment::assayNames(sce))
  expect_equal(nrow(SummarizedExperiment::colData(sce)), N_CELLS)
})

test_that("LoadH5SCE counts values match LoadH5Seurat counts", {
  skip_if_no_sce()
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  sce  <- LoadH5SCE(path, verbose = FALSE)
  seur <- LoadH5Seurat(path, verbose = FALSE)

  expect_equal(
    as.matrix(SummarizedExperiment::assay(sce, "counts")),
    as.matrix(SeuratObject::GetAssayData(seur, layer = "counts"))
  )
})

# ---------------------------------------------------------------------------
# cells= subset
# ---------------------------------------------------------------------------

test_that("LoadH5SCE cells= character loads only requested cells", {
  skip_if_no_sce()
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  all_cells  <- colnames(make_test_seurat())
  keep       <- all_cells[c(1L, 5L, 10L, 15L, 20L)]
  sce        <- LoadH5SCE(path, cells = keep, verbose = FALSE)

  expect_equal(ncol(sce), length(keep))
  expect_equal(colnames(sce), keep)
  expect_equal(rownames(SummarizedExperiment::colData(sce)), keep)
})

test_that("LoadH5SCE cells= subset matches manually-subset full load", {
  skip_if_no_sce()
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  all_cells  <- colnames(make_test_seurat())
  keep       <- all_cells[c(3L, 7L, 11L)]

  sce_full <- LoadH5SCE(path, verbose = FALSE)
  sce_sub  <- LoadH5SCE(path, cells = keep, verbose = FALSE)

  expect_equal(
    as.matrix(SummarizedExperiment::assay(sce_sub, "counts")),
    as.matrix(SummarizedExperiment::assay(sce_full, "counts")[, keep])
  )
})

test_that("LoadH5SCE rejects unknown cell names", {
  skip_if_no_sce()
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  expect_error(
    LoadH5SCE(path, cells = c("cell001", "GHOST"), verbose = FALSE),
    regexp = "GHOST"
  )
})

# ---------------------------------------------------------------------------
# as.delayed = TRUE
# ---------------------------------------------------------------------------

test_that("LoadH5SCE(as.delayed = TRUE) returns a DelayedArray-backed counts", {
  skip_if_no_sce()
  testthat::skip_if_not_installed("HDF5Array")
  testthat::skip_if_not_installed("DelayedArray")
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  sce <- LoadH5SCE(path, as.delayed = TRUE, verbose = FALSE)
  m   <- SummarizedExperiment::assay(sce, "counts")

  expect_true(is(m, "DelayedArray"))
  expect_equal(dim(m), c(N_FEATURES, N_CELLS))

  sce_eager <- LoadH5SCE(path, verbose = FALSE)
  expect_equal(
    as.matrix(m),
    as.matrix(SummarizedExperiment::assay(sce_eager, "counts"))
  )
})

test_that("LoadH5SCE(as.delayed = TRUE) supports cells= via lazy subset", {
  skip_if_no_sce()
  testthat::skip_if_not_installed("HDF5Array")
  testthat::skip_if_not_installed("DelayedArray")
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  all_cells <- colnames(make_test_seurat())
  keep      <- all_cells[c(2L, 4L, 6L)]

  sce_full <- LoadH5SCE(path, verbose = FALSE)
  sce_lazy <- LoadH5SCE(path, cells = keep, as.delayed = TRUE, verbose = FALSE)

  expect_true(is(SummarizedExperiment::assay(sce_lazy, "counts"), "DelayedArray"))
  expect_equal(ncol(sce_lazy), length(keep))
  expect_equal(
    as.matrix(SummarizedExperiment::assay(sce_lazy, "counts")),
    as.matrix(SummarizedExperiment::assay(sce_full, "counts")[, keep])
  )
})
