#' @include zzz.R
#' @include h5Seurat.R
#' @include ReadH5.R
#'
NULL

#' Load an h5Seurat file as a \code{SingleCellExperiment}
#'
#' Reads counts and (optionally) cell-level metadata from an h5Seurat file and
#' returns a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.  This is
#' an initial implementation and intentionally narrow in scope: only the
#' \code{counts} slot of each requested assay and the cell metadata are
#' loaded.  All other components (\code{reductions}, \code{graphs},
#' \code{neighbors}, \code{images}, \code{commands}, \code{misc}, \code{tools})
#' are accepted as parameters but must be set to \code{FALSE}; passing any
#' other value raises an error so that callers do not silently get an
#' incomplete object.
#'
#' @section DelayedArray (\code{as.delayed}):
#' When \code{as.delayed = TRUE} the counts matrix is wrapped as a
#' \code{\link[DelayedArray]{DelayedArray}} backed by
#' \code{\link[HDF5Array]{H5SparseMatrixSeed}} and never realised in memory.
#' Cell subsetting via \code{cells} is then a lazy
#' \code{DelayedSubset}; rows and columns are only materialised when the
#' caller pulls values out of the matrix.
#'
#' Implementation notes (research outcome):
#' \itemize{
#'   \item h5seurat stores counts as a CSC group with datasets \code{data},
#'     \code{indices}, \code{indptr} and a \code{dims} attribute.
#'     \code{HDF5Array} cannot auto-detect this layout (it expects 10x's
#'     \code{/matrix/shape} or AnnData's \code{h5sparse_format} markers), so
#'     this loader passes \code{dim} and \code{sparse.layout = "csc"}
#'     explicitly to \code{H5SparseMatrixSeed()}.
#'   \item The \code{hdf5r}/\code{rhdf5} file-lock conflict that killed
#'     \href{https://github.com/mojaveazure/seurat-disk/pull/37}{PR #37} arose
#'     because both libraries were holding the file open for writes.  For
#'     read-only loading we close the \code{hdf5r} connection before handing
#'     the path to \code{HDF5Array}, which avoids the conflict.
#'   \item \code{HDF5Array} is an optional dependency.  When
#'     \code{as.delayed = TRUE} but the package is not installed, an error is
#'     raised pointing the user at Bioconductor.
#' }
#'
#' @param file Path to an h5Seurat file
#' @param assay Name of the assay to load, or \code{NULL} for the default
#'   assay.  Only a single assay is supported in this initial implementation
#' @param cells \code{NULL} (load all cells), or a logical, integer, or
#'   character vector identifying the cells to load.  See
#'   \code{\link{LoadH5Seurat}} for the resolution rules
#' @param reductions,graphs,neighbors,images,commands,misc,tools Must be
#'   \code{FALSE}; reserved for future expansion
#' @param meta.data Whether to load cell-level metadata into \code{colData}
#' @param as.delayed If \code{TRUE}, keep counts on disk as a
#'   \code{\link[DelayedArray]{DelayedArray}}; requires the
#'   \code{HDF5Array} package
#' @param verbose Show progress messages
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'
#' @export
#'
LoadH5SCE <- function(
  file,
  assay       = NULL,
  cells       = NULL,
  reductions  = FALSE,
  graphs      = FALSE,
  neighbors   = FALSE,
  images      = FALSE,
  meta.data   = TRUE,
  commands    = FALSE,
  misc        = FALSE,
  tools       = FALSE,
  as.delayed  = FALSE,
  verbose     = TRUE
) {
  if (!requireNamespace('SingleCellExperiment', quietly = TRUE)) {
    stop(
      "LoadH5SCE requires the 'SingleCellExperiment' Bioconductor package; ",
      "install with BiocManager::install('SingleCellExperiment')",
      call. = FALSE
    )
  }
  if (!is.character(file) || length(file) != 1L) {
    stop("'file' must be a single path to an h5Seurat file", call. = FALSE)
  }
  unsupported <- list(
    reductions = reductions,
    graphs     = graphs,
    neighbors  = neighbors,
    images     = images,
    commands   = commands,
    misc       = misc,
    tools      = tools
  )
  bad <- vapply(unsupported, FUN = function(v) !isFALSE(v), FUN.VALUE = logical(1L))
  if (any(bad)) {
    stop(
      "LoadH5SCE only supports counts and metadata in this initial ",
      "implementation; the following arguments must be FALSE: ",
      paste(names(unsupported)[bad], collapse = ", "),
      call. = FALSE
    )
  }

  hfile <- h5Seurat$new(filename = file, mode = 'r')
  on.exit(expr = if (hfile$is_valid) hfile$close_all(), add = TRUE)

  index           <- hfile$index()
  index.assays    <- setdiff(names(index), c('global', 'no.assay'))
  if (is.null(assay)) {
    assay <- DefaultAssay(object = hfile)
  } else if (!is.character(assay) || length(assay) != 1L) {
    stop("'assay' must be NULL or a single assay name", call. = FALSE)
  }
  if (!assay %in% index.assays) {
    stop("Assay not found in file: ", assay, call. = FALSE)
  }

  all.cells       <- Cells(x = hfile)
  cell_idx        <- resolve_cells(cells = cells, all_cells = all.cells)
  selected.cells  <- if (is.null(cell_idx)) all.cells else all.cells[cell_idx]

  grp <- hfile[['assays']][[assay]]
  if (!grp$exists(name = 'counts')) {
    stop("Assay '", assay, "' has no 'counts' slot to load", call. = FALSE)
  }
  counts.obj <- grp[['counts']]
  if (as.delayed && !inherits(counts.obj, 'H5Group')) {
    stop("as.delayed = TRUE currently only supports sparse (CSC) counts; ",
         "assay '", assay, "' stores 'counts' as a dense dataset", call. = FALSE)
  }
  features <- FixFeatures(features = grp[['features']][])

  # Pull what we need from the file before any HDF5Array call: HDF5Array
  # re-opens the file via rhdf5 and cannot share file locks with the live
  # hdf5r handle.
  counts.dims <- NULL
  counts.path <- NULL
  counts.mat  <- NULL
  if (as.delayed) {
    if (!counts.obj$attr_exists(attr_name = 'dims')) {
      stop("Counts group has no 'dims' attribute; cannot wrap as DelayedArray",
           call. = FALSE)
    }
    counts.dims <- hdf5r::h5attr(x = counts.obj, which = 'dims')
    counts.path <- counts.obj$get_obj_name()
  } else {
    if (verbose) {
      message("Reading counts for assay '", assay, "'")
    }
    counts.mat <- read_assay_matrix(x = counts.obj, col_idx = cell_idx)
    rownames(counts.mat) <- features
    colnames(counts.mat) <- selected.cells
  }

  # Read cell metadata before closing the handle.
  col.data <- NULL
  if (isTRUE(meta.data) && hfile$exists(name = 'meta.data')) {
    if (verbose) {
      message("Reading cell metadata")
    }
    md <- as.data.frame(x = hfile[['meta.data']], row.names = all.cells)
    if (!is.null(cell_idx)) {
      md <- md[cell_idx, , drop = FALSE]
    }
    col.data <- S4Vectors::DataFrame(md, check.names = FALSE)
  }

  if (as.delayed) {
    # Materialise everything we need from the file and release the hdf5r
    # handle so HDF5Array can reopen with rhdf5.
    abs.path <- normalizePath(file, mustWork = TRUE)
    hfile$close_all()
    if (!requireNamespace('HDF5Array', quietly = TRUE)) {
      stop(
        "as.delayed = TRUE requires the 'HDF5Array' Bioconductor package; ",
        "install with BiocManager::install('HDF5Array')",
        call. = FALSE
      )
    }
    if (!requireNamespace('DelayedArray', quietly = TRUE)) {
      stop(
        "as.delayed = TRUE requires the 'DelayedArray' Bioconductor package",
        call. = FALSE
      )
    }
    seed <- HDF5Array::H5SparseMatrixSeed(
      filepath      = abs.path,
      group         = counts.path,
      dim           = as.integer(counts.dims),
      sparse.layout = "csc"
    )
    counts.mat <- DelayedArray::DelayedArray(seed = seed)
    dimnames(counts.mat) <- list(features, all.cells)
    if (!is.null(cell_idx)) {
      counts.mat <- counts.mat[, cell_idx, drop = FALSE]
    }
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts.mat)
  )
  if (!is.null(col.data)) {
    SummarizedExperiment::colData(sce) <- col.data
  }
  return(sce)
}
