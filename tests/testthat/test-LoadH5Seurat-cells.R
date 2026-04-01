# Tests for LoadH5Seurat cells= parameter
#
# Organisation
#   Baseline (§1) – roundtrip correctness BEFORE cells= exists; these must
#                   pass against the unmodified package too.
#   cells= feature  – verify HDF5-level subsetting once implemented.
#   resolve_cells() – unit tests for the internal helper; no HDF5 needed.

# ---------------------------------------------------------------------------
# §1  Baseline roundtrip tests
# ---------------------------------------------------------------------------

test_that("baseline: full roundtrip gives correct cell and feature counts", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  obj <- LoadH5Seurat(path, verbose = FALSE)

  expect_equal(ncol(obj), N_CELLS)
  expect_equal(nrow(obj), N_FEATURES)
})

test_that("baseline: full roundtrip gives correct PCA embedding shape", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  obj <- LoadH5Seurat(path, verbose = FALSE)
  emb <- Seurat::Embeddings(obj, "pca")

  expect_equal(nrow(emb), N_CELLS)
  expect_equal(ncol(emb), N_PCS)
})

test_that("baseline: full roundtrip preserves cell names", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  expected_cells <- colnames(make_test_seurat())
  obj <- LoadH5Seurat(path, verbose = FALSE)

  expect_equal(colnames(obj), expected_cells)
})

test_that("baseline: full roundtrip metadata row count matches ncol", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  obj <- LoadH5Seurat(path, verbose = FALSE)

  expect_equal(nrow(obj[[]]), N_CELLS)
})

test_that("baseline: manual subset of loaded object gives correct shapes", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  obj_full  <- LoadH5Seurat(path, verbose = FALSE)
  n_keep    <- 5L
  keep_cells <- colnames(obj_full)[seq_len(n_keep)]
  obj_sub   <- obj_full[, keep_cells]

  expect_equal(ncol(obj_sub), n_keep)
  expect_equal(colnames(obj_sub), keep_cells)
  expect_equal(nrow(Seurat::Embeddings(obj_sub, "pca")), n_keep)
  expect_equal(nrow(obj_sub[[]]), n_keep)
})

# ---------------------------------------------------------------------------
# §2  resolve_cells() unit tests  (pure R, no HDF5)
# ---------------------------------------------------------------------------

test_that("resolve_cells: NULL returns NULL", {
  expect_null(SeuratDisk:::resolve_cells(NULL, letters))
})

test_that("resolve_cells: logical vector returns which()", {
  all_cells <- paste0("c", 1:10)
  mask <- c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
  expect_equal(SeuratDisk:::resolve_cells(mask, all_cells), which(mask))
})

test_that("resolve_cells: logical vector wrong length errors", {
  expect_error(
    SeuratDisk:::resolve_cells(c(TRUE, FALSE), paste0("c", 1:10)),
    regexp = "length"
  )
})

test_that("resolve_cells: integer vector passes through as integer", {
  all_cells <- paste0("c", 1:10)
  expect_equal(SeuratDisk:::resolve_cells(c(1, 3, 5), all_cells), c(1L, 3L, 5L))
})

test_that("resolve_cells: integer out of range errors", {
  expect_error(
    SeuratDisk:::resolve_cells(c(1L, 99L), paste0("c", 1:10)),
    regexp = "out of range"
  )
})

test_that("resolve_cells: character vector returns correct indices", {
  all_cells <- paste0("cell", 1:20)
  query     <- c("cell3", "cell7", "cell15")
  expect_equal(SeuratDisk:::resolve_cells(query, all_cells), c(3L, 7L, 15L))
})

test_that("resolve_cells: unknown barcode errors with its name", {
  all_cells <- paste0("cell", 1:5)
  err <- expect_error(
    SeuratDisk:::resolve_cells(c("cell1", "GHOST"), all_cells)
  )
  expect_match(conditionMessage(err), "GHOST")
})

# ---------------------------------------------------------------------------
# §3  cells= feature tests  (require cells= implementation)
# ---------------------------------------------------------------------------

test_that("cells= character: loads only requested cells", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  all_cells  <- colnames(make_test_seurat())
  keep_cells <- all_cells[c(1L, 5L, 10L, 15L, 20L)]

  obj <- LoadH5Seurat(path, cells = keep_cells, verbose = FALSE)

  expect_equal(ncol(obj), length(keep_cells))
  expect_equal(sort(colnames(obj)), sort(keep_cells))
})

test_that("cells= integer: loads only requested cells", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  idx        <- c(2L, 4L, 6L, 8L, 10L)
  all_cells  <- colnames(make_test_seurat())
  keep_cells <- all_cells[idx]

  obj <- LoadH5Seurat(path, cells = idx, verbose = FALSE)

  expect_equal(ncol(obj), length(idx))
  expect_equal(sort(colnames(obj)), sort(keep_cells))
})

test_that("cells= logical: loads only TRUE cells", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  mask       <- rep(c(TRUE, FALSE), times = N_CELLS / 2)
  all_cells  <- colnames(make_test_seurat())
  keep_cells <- all_cells[mask]

  obj <- LoadH5Seurat(path, cells = mask, verbose = FALSE)

  expect_equal(ncol(obj), sum(mask))
  expect_equal(sort(colnames(obj)), sort(keep_cells))
})

test_that("cells= NULL: behaviour unchanged (all cells loaded)", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  obj <- LoadH5Seurat(path, cells = NULL, verbose = FALSE)

  expect_equal(ncol(obj), N_CELLS)
})

test_that("cells=: PCA embedding rows match subset", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  keep_cells <- colnames(make_test_seurat())[1:5]
  obj <- LoadH5Seurat(path, cells = keep_cells, verbose = FALSE)
  emb <- Seurat::Embeddings(obj, "pca")

  expect_equal(nrow(emb), length(keep_cells))
  expect_equal(rownames(emb), keep_cells)
})

test_that("cells=: metadata rows match subset", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  keep_cells <- colnames(make_test_seurat())[1:7]
  obj <- LoadH5Seurat(path, cells = keep_cells, verbose = FALSE)

  expect_equal(nrow(obj[[]]), length(keep_cells))
  expect_equal(rownames(obj[[]]), keep_cells)
})

test_that("cells=: counts matrix columns match subset", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  keep_cells <- colnames(make_test_seurat())[c(3L, 7L, 11L)]
  obj    <- LoadH5Seurat(path, cells = keep_cells, verbose = FALSE)
  obj_full <- LoadH5Seurat(path, verbose = FALSE)

  expect_equal(
    as.matrix(SeuratObject::GetAssayData(obj, layer = "counts")),
    as.matrix(SeuratObject::GetAssayData(obj_full, layer = "counts")[, keep_cells])
  )
})

test_that("cells=: scale.data columns match subset", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  keep_cells <- colnames(make_test_seurat())[c(2L, 8L, 14L)]
  obj      <- LoadH5Seurat(path, cells = keep_cells, verbose = FALSE)
  obj_full <- LoadH5Seurat(path, verbose = FALSE)

  expect_equal(
    SeuratObject::GetAssayData(obj, layer = "scale.data"),
    SeuratObject::GetAssayData(obj_full, layer = "scale.data")[, keep_cells]
  )
})

test_that("cells=: PCA values match manual subset", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  keep_cells <- colnames(make_test_seurat())[c(1L, 3L, 5L)]
  obj      <- LoadH5Seurat(path, cells = keep_cells, verbose = FALSE)
  obj_full <- LoadH5Seurat(path, verbose = FALSE)

  expect_equal(
    Seurat::Embeddings(obj, "pca"),
    Seurat::Embeddings(obj_full, "pca")[keep_cells, ]
  )
})

test_that("cells=: unknown barcode errors informatively", {
  path <- make_test_h5seurat()
  withr::defer(unlink(path))

  expect_error(
    LoadH5Seurat(path, cells = c("cell001", "GHOST_CELL"), verbose = FALSE),
    regexp = "GHOST_CELL"
  )
})
