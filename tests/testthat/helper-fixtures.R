# Force v3-style Assay objects so SaveH5Seurat (written for Seurat ≤4) works.
# Seurat 5 defaults to Assay5 which uses layers instead of slots; SeuratDisk
# has not been updated for Assay5 yet.
options(Seurat.object.assay.version = "v3")

# Shared constants -----------------------------------------------------------

N_CELLS    <- 20L
N_FEATURES <- 50L
N_PCS      <- 5L

# Fixture builders -----------------------------------------------------------

#' Create a minimal in-memory Seurat object for testing.
#'
#' Includes: sparse counts, normalised data, dense scale.data (all features),
#' a manually-constructed PCA reduction, and two metadata columns.
make_test_seurat <- function() {
  set.seed(42L)

  counts_mat <- abs(matrix(
    data = rnorm(n = N_FEATURES * N_CELLS, mean = 3, sd = 1),
    nrow = N_FEATURES,
    ncol = N_CELLS
  ))
  counts_mat <- round(counts_mat)
  rownames(counts_mat) <- paste0("gene", seq_len(N_FEATURES))
  colnames(counts_mat) <- paste0("cell", formatC(seq_len(N_CELLS), width = 3L, flag = "0"))
  counts <- Matrix::Matrix(counts_mat, sparse = TRUE)

  obj <- SeuratObject::CreateSeuratObject(
    counts      = counts,
    min.cells   = 0L,
    min.features = 0L
  )
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, features = rownames(obj), verbose = FALSE)

  # Build a PCA without running irlba so tests stay fast and dependency-light
  pca_emb <- matrix(rnorm(N_CELLS * N_PCS), nrow = N_CELLS, ncol = N_PCS)
  rownames(pca_emb) <- colnames(counts)
  colnames(pca_emb) <- paste0("PC_", seq_len(N_PCS))
  obj[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = pca_emb,
    key        = "PC_",
    assay      = "RNA"
  )

  obj$condition <- factor(sample(c("ctrl", "stim"), N_CELLS, replace = TRUE))
  obj$score     <- runif(N_CELLS)

  return(obj)
}

#' Save make_test_seurat() to a temp .h5seurat file and return the path.
#' The caller is responsible for cleanup (use withr::defer(unlink(path))).
make_test_h5seurat <- function() {
  path <- tempfile(fileext = ".h5seurat")
  SeuratDisk::SaveH5Seurat(
    object   = make_test_seurat(),
    filename = path,
    overwrite = TRUE,
    verbose  = FALSE
  )
  return(path)
}
