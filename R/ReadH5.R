#' @include TestObject.R TestH5.R
#' @importFrom hdf5r h5attr
#'
NULL

#' Load data from an HDF5 File
#'
#' HDF5 allows storing data in an arbitrary fashion, which makes reading data
#' into memory a hassle. The methods here serve as convenience functions for
#' reading data stored in a certain format back into a certain R object. For
#' details regarding how data should be stored on disk, please see the
#' \href{../doc/h5Seurat-spec.html}{h5Seurat file specification}.
#'
#' @param x An HDF5 dataset or group
#' @param ... Arguments passed to other methods
#'
#' @name ReadH5
#' @rdname ReadH5
#'
#' @md
#'
NULL

#' @return \code{as.array}: returns an \code{\link[base]{array}} with the data
#' from the HDF5 dataset
#'
#' @aliases as.array
#'
#' @rdname ReadH5
#' @method as.array H5D
#' @export
#'
as.array.H5D <- function(x, ...) {
  return(as.array(x = x$read(), ...))
}

#' @inheritParams base::as.data.frame
#'
#' @return \code{as.data.frame}: returns a \code{\link[base]{data.frame}} with
#' the data from the HDF5 dataset or group
#'
#' @aliases as.data.frame
#'
#' @rdname ReadH5
#' @method as.data.frame H5D
#' @export
#'
as.data.frame.H5D <- function(x, row.names = NULL, optional = FALSE, ...) {
  df <- x[]
  if (!is.null(x = row.names)) {
    row.names(x = df) <- row.names
  }
  if (!optional) {
    colnames(x = df) <- make.names(names = x$get_type()$get_cpd_labels())
  }
  if (x$attr_exists(attr_name = 'logicals')) {
    bool.cols <- intersect(
      x = h5attr(x = x, which = 'logicals'),
      y = colnames(x = df)
    )
    for (i in bool.cols) {
      bools <- df[, i, drop = TRUE]
      bools[bools == 2] <- NA_integer_
      bools <- as.logical(x = bools)
      df[[i]] <- bools
    }
  }
  return(df)
}

#' @rdname ReadH5
#' @method as.data.frame H5Group
#' @export
#'
as.data.frame.H5Group <- function(x, row.names = NULL, optional = FALSE, ...) {
  df.names <- names(x = x)
  if (x$attr_exists(attr_name = 'colnames')) {
    df.order <- h5attr(x = x, which = 'colnames')
    missing.cols <- setdiff(x = df.order, y = df.names)
    if (length(x = missing.cols)) {
      if (length(x = missing.cols) == length(x = df.order)) {
        warning(
          "None of the columns specified by 'colnames' are present",
          call. = FALSE,
          immediate. = TRUE
        )
        df.order <- df.names
      } else {
        warning(
          "The following columns specified by 'colnames' are missing: ",
          paste(missing.cols, collapse = ', '),
          call. = FALSE,
          immediate. = TRUE
        )
        df.order <- setdiff(x = df.order, y = missing.cols)
      }
    }
    df.names <- c(df.order, df.names[!df.names %in% df.order])
  }
  df <- vector(mode = 'list', length = length(x = df.names))
  names(x = df) <- df.names
  bool.cols <- if (x$attr_exists(attr_name = 'logicals')) {
    intersect(x = h5attr(x = x, which = 'logicals'), y = colnames(x = df))
  } else {
    character(length = 0L)
  }
  for (i in df.names) {
    if (inherits(x = x[[i]], what = 'H5D')) {
      df[[i]] <- if (i %in% bool.cols) {
        as.logical(x = x[[i]])
      } else {
        x[[i]][]
      }
    } else if (inherits(x = x[[i]], what = 'H5Group')) {
      # Call the internal helper directly: relying on as.factor() S4 dispatch
      # here is fragile because as.factor is implicitly promoted to a generic
      # and other packages on the search path (BiocGenerics / S4Vectors stack
      # pulled in by SingleCellExperiment, scater, etc.) can shadow our
      # generic, dropping the H5Group method and falling through to base.
      df[[i]] <- h5group_to_factor(x = x[[i]])
    } else {
      stop("Unknown dataset type for ", i, call. = FALSE)
    }
  }
  if (AttrExists(x = x, name = '_index')) {
    rnames <- h5attr(x = x, which = '_index')
    if (is.null(x = row.names)) {
      row.names <- x[[rnames]][]
    }
    df[[rnames]] <- NULL
    df.names <- setdiff(x = df.names, rnames)
  }
  df <- as.data.frame(x = df, row.names = row.names, optional = optional, ...)
  names(x = df) <- df.names
  return(df)
}

#' @importFrom stats na.omit
#' @importFrom methods setMethod
#'
#' @return \code{as.factor}: returns a \code{\link[base]{factor}} with the data
#' from the HDF5 group
#'
#' @aliases as.factor
#'
#' @rdname ReadH5
#' @export
#'
#' Read an H5Group encoding a factor (datasets \code{levels} + \code{values})
#' as an R \code{\link[base]{factor}}.
#'
#' Internal helper used in place of \code{as.factor(<H5Group>)} S4 dispatch,
#' which is fragile across package load orders (see ReadH5.R).
#'
#' @keywords internal
h5group_to_factor <- function(x) {
  if (!x$exists(name = 'levels') || !x$exists(name = 'values')) {
    stop("Missing required datasets 'levels' and 'values'", call. = FALSE)
  }
  if (!IsDType(x = x[['levels']], dtype = 'H5T_STRING') || length(x = x[['levels']]$dims) != 1) {
    stop("'levels' must be a one-dimensional string dataset", call. = FALSE)
  }
  if (!IsDType(x = x[['values']], dtype = 'H5T_INTEGER') || length(x = x[['values']]$dims) != 1) {
    stop("'values' must be a one-dimensional integer dataset", call. = FALSE)
  }
  if (!x[['levels']]$dims) {
    return(factor())
  }
  values <- x[['values']][]
  levels <- x[['levels']][]
  if (length(x = unique(x = na.omit(object = values))) > length(x = levels)) {
    stop("Too many values for levels provided", call. = FALSE)
  }
  factor(x = levels[values], levels = levels)
}

setMethod(
  f = 'as.factor',
  signature = c('x' = 'H5Group'),
  definition = function(x) h5group_to_factor(x = x)
)

#' @importFrom withr with_package
#'
#' @return \code{as.list}: returns a \code{\link[base]{list}} with the data from
#' the HDF5 group
#'
#' @aliases as.list
#'
#' @rdname ReadH5
#' @export
#'
setMethod(
  f = 'as.list',
  signature = c('x' = 'H5Group'),
  definition = function(x, which = NULL, ...) {
    list.names <- which %||% names(x = x)
    is.vec.name <- grepl(
      pattern = "__names__",
      x = list.names,
      ignore.case = FALSE,
      perl = FALSE,
      fixed = TRUE
    )
    vec.name.obj <- list.names[is.vec.name]
    list.names <- list.names[!is.vec.name]
    if (x$attr_exists(attr_name = 'names')) {
      list.order <- h5attr(x = x, which = 'names')
      missing.names <- setdiff(x = list.order, y = list.names)
      if (length(x = missing.names)) {
        if (length(x = missing.names) == length(x = list.order)) {
          warning(
            "None of the named entires specified by 'names' are present",
            call. = FALSE,
            immediate. = TRUE
          )
          list.order <- list.names
        } else {
          warning(
            "The following named entries specified by 'names' are missing: ",
            paste(missing.names, collapse = ', '),
            call. = FALSE,
            immediate. = TRUE
          )
          list.order <- setdiff(x = list.order, y = missing.names)
        }
      }
      list.names <- c(list.order, list.names[!list.names %in% list.order])
    }
    data <- vector(mode = 'list', length = length(x = list.names))
    names(x = data) <- list.names
    for (i in list.names) {
      i_name <- paste0(i, "__names__")
      if (inherits(x = x[[i]], what = 'H5D')) {
        data[[i]] <- if (IsDataFrame(x = x[[i]])) {
          as.data.frame(x = x[[i]], ...)
        } else if (IsMatrix(x = x[[i]])) {
          as.matrix(x = x[[i]], ...)
        } else if (IsLogical(x = x[[i]])) {
          v <- as.logical(x = x[[i]], ...)
          if (i_name %in% vec.name.obj) {
            names(x = v) <- x[[i_name]]$read()
          }
          v
        } else if (!x[[i]]$dims) {
          NULL
        } else {
          v <- x[[i]]$read()
          if (i_name %in% vec.name.obj) {
            names(x = v) <- x[[i_name]]$read()
          }
          v
        }
      } else {
        data[[i]] <- if (IsFactor(x = x[[i]])) {
          as.factor(x = x[[i]])
        } else if (IsDataFrame(x = x[[i]])) {
          as.data.frame(x = x[[i]], ...)
        } else if (IsMatrix(x = x[[i]])) {
          as.matrix(x = x[[i]], ...)
        } else {
          as.list(x = x[[i]], ...)
        }
      }
    }
    if (x$attr_exists(attr_name = 's3class')) {
      data <- structure(.Data = data, class = h5attr(x = x, which = 's3class'))
    } else if (x$attr_exists(attr_name = 's4class')) {
      attr(x = data, which = 'classDef') <- h5attr(x = x, which = 's4class')
      data <- SeuratObject::ListToS4(x = data)
      # if (grepl(pattern = ':', x = class)) {
      #   classdef <- unlist(x = strsplit(x = class, split = ':'))
      #   classpkg <- classdef[1]
      #   class <- classdef[2]
      #   try(
      #     expr = class <- with_package(
      #       package = classpkg,
      #       code = getClass(Class = class)
      #     ),
      #     silent = TRUE
      #   )
      # }
      # try(
      #   expr = data <- do.call(what = 'new', args = c('Class' = class, data)),
      #   silent = TRUE
      # )
    }
    return(data)
  }
)

#' @return \code{as.logical}: returns a \code{\link[base]{logical}} with the
#' data from the HDF5 dataset
#'
#' @aliases as.logical
#'
#' @rdname ReadH5
#' @method as.logical H5D
#' @export
#'
as.logical.H5D <- function(x, ...) {
  bool <- x$read()
  bool[which(x = bool == 2)] <- NA_integer_
  return(as.logical(x = bool, ...))
}

#' @param transpose Transpose the data upon reading it in, used when writing
#' data in row-major order (eg. from C or Python)
#'
#' @return \code{as.matrix}, \code{H5D} method: returns a
#' \code{\link[base]{matrix}} with the data from the HDF5 dataset
#'
#' @aliases as.matrix
#'
#' @rdname ReadH5
#' @method as.matrix H5D
#' @export
#'
as.matrix.H5D <- function(x, transpose = FALSE, ...) {
  obj <- x$read()
  if (transpose) {
    obj <- t(x = obj)
  }
  return(as.matrix(x = obj))
}

#' @rdname ReadH5
#' @method as.matrix H5Group
#' @export
#'
as.matrix.H5Group <- function(x, ...) {
  return(as.sparse(x = x, ...))
}

#' @importFrom Matrix Matrix
#' @importFrom Seurat as.sparse
#' @importFrom utils setTxtProgressBar
#'
#' @return \code{as.sparse}, \code{H5D} method: returns a sparse matrix with the
#' data from the HDF5 dataset
#'
#' @rdname ReadH5
#' @method as.sparse H5D
#' @export
#'
as.sparse.H5D <- function(x, verbose = TRUE, ...) {
  xdims <- x$dims
  MARGIN <- GetMargin(dims = xdims)
  chunk.points <- ChunkPoints(
    dsize = xdims[MARGIN],
    csize = x$chunk_dims[MARGIN]
  )
  mat <- Matrix(data = 0, nrow = xdims[1], ncol = xdims[2], sparse = TRUE)
  dims <- vector(mode = 'list', length = 2L)
  dims[[-MARGIN]] <- seq_len(length.out = xdims[-MARGIN])
  if (isTRUE(x = verbose)) {
    pb <- PB()
  }
  for (i in seq_len(length.out = nrow(x = chunk.points))) {
    dims[[MARGIN]] <- seq.default(
      from = chunk.points[i, 'start'],
      to = chunk.points[i, 'end']
    )
    mat[dims[[1]], dims[[2]]] <- x$read(args = dims, drop = FALSE)
    if (isTRUE(x = verbose)) {
      setTxtProgressBar(pb = pb, value = i / nrow(x = chunk.points))
    }
  }
  if (isTRUE(x = verbose)) {
    close(con = pb)
  }
  return(as.sparse(x = mat))
}

#' @importFrom Seurat as.sparse
#' @importFrom Matrix sparseMatrix
#'
#' @return \code{as.sparse}, \code{as.matrix}, \code{H5Group} method: returns a
#' \code{\link[Matrix]{sparseMatrix}} with the data from the HDF5 group
#'
#' @aliases as.sparse
#'
#' @rdname ReadH5
#' @method as.sparse H5Group
#' @export
#'
as.sparse.H5Group <- function(x, ...) {
  for (i in c('data', 'indices', 'indptr')) {
    if (!x$exists(name = i) || !inherits(x = x[[i]], what = 'H5D')) {
      stop("Missing dataset ", i, call. = FALSE)
    } else if (length(x = x[[i]]$dims) != 1) {
      stop("Dataset ", i, " is not one-dimensional", call. = FALSE)
    }
    if (i == 'data' && !IsDType(x = x[[i]], dtype = c('H5T_FLOAT', 'H5T_INTEGER'))) {
      stop("'data' must be integer or floating-point values", call. = FALSE)
    } else if (i != 'data' && !IsDType(x = x[[i]], dtype = 'H5T_INTEGER')) {
      stop("'", i, "' must be integer values", call. = FALSE)
    }
  }
  if (x$attr_exists(attr_name = 'dims')) {
    return(sparseMatrix(
      i = x[['indices']][] + 1,
      p = x[['indptr']][],
      x = x[['data']][],
      dims = h5attr(x = x, which = 'dims')
    ))
  }
  # Because AnnData likes changing stuff
  adata.shape <- vapply(
    X = c('h5sparse_shape', 'shape'),
    FUN = x$attr_exists,
    FUN.VALUE = logical(length = 1L)
  )
  if (any(adata.shape)) {
    return(sparseMatrix(
      i = x[['indices']][] + 1,
      p = x[['indptr']][],
      x = x[['data']][],
      dims = rev(x = h5attr(
        x = x,
        which = names(x = which(x = adata.shape))[1]
      ))
    ))
  }
  return(sparseMatrix(
    i = x[['indices']][] + 1,
    p = x[['indptr']][],
    x = x[['data']][]
  ))
}

#' Read specific columns from a CSC sparse matrix stored in an HDF5 group
#'
#' Performs HDF5-level column slicing on a group containing \code{data},
#' \code{indices}, and \code{indptr} datasets (standard CSC layout).  Only
#' the entries for the requested columns are read from disk; the full value
#' and index arrays are never realised in memory.
#'
#' @param grp An \code{H5Group} with CSC datasets and a \code{dims}
#'   (or \code{h5sparse_shape} / \code{shape}) attribute.
#' @param col_idx Integer vector of 1-based column indices to extract.
#'
#' @return A \code{\link[Matrix]{dgCMatrix}} with \code{length(col_idx)}
#'   columns.
#'
#' @importFrom Matrix sparseMatrix
#'
#' @keywords internal
#'
read_sparse_cols <- function(grp, col_idx) {
  indptr_ds  <- grp[['indptr']]
  indices_ds <- grp[['indices']]
  data_ds    <- grp[['data']]

  # Resolve matrix dimensions from whichever attribute is present
  if (grp$attr_exists(attr_name = 'dims')) {
    dims <- h5attr(x = grp, which = 'dims')
  } else if (grp$attr_exists(attr_name = 'h5sparse_shape')) {
    dims <- rev(h5attr(x = grp, which = 'h5sparse_shape'))
  } else if (grp$attr_exists(attr_name = 'shape')) {
    dims <- rev(h5attr(x = grp, which = 'shape'))
  } else {
    stop("Cannot determine matrix dimensions: no dims/shape attribute", call. = FALSE)
  }
  nrow_out   <- dims[1L]
  ncol_out   <- length(col_idx)

  # Read the full indptr array; it is small (ncols * 4 bytes) and needed for
  # every column's byte-range calculation.
  indptr <- indptr_ds$read()

  # For column c (1-based R):
  #   data start (1-based) = indptr[c]   + 1   (indptr is 0-based in file)
  #   data end   (1-based) = indptr[c+1]        (0-based exclusive = 1-based inclusive)
  starts <- indptr[col_idx]        + 1L   # 1-based start positions
  ends   <- indptr[col_idx + 1L]          # 1-based end positions (inclusive)

  row_idx_list <- vector('list', ncol_out)
  val_list     <- vector('list', ncol_out)
  new_indptr   <- integer(ncol_out + 1L)
  new_indptr[1L] <- 0L

  for (j in seq_len(ncol_out)) {
    s <- starts[j]
    e <- ends[j]
    if (e >= s) {
      rng <- s:e
      row_idx_list[[j]] <- indices_ds$read(args = list(rng)) + 1L
      val_list[[j]]     <- data_ds$read(args = list(rng))
      new_indptr[j + 1L] <- new_indptr[j] + (e - s + 1L)
    } else {
      new_indptr[j + 1L] <- new_indptr[j]
    }
  }

  all_row_idx <- unlist(row_idx_list, use.names = FALSE)
  all_vals    <- unlist(val_list,     use.names = FALSE)

  sparseMatrix(
    i    = if (length(all_row_idx)) all_row_idx else integer(0L),
    p    = new_indptr,
    x    = if (length(all_vals))    all_vals    else numeric(0L),
    dims = c(nrow_out, ncol_out)
  )
}

#' Read a matrix slot from an HDF5 assay dataset, optionally subsetting columns
#'
#' Dispatches on the HDF5 object type:
#' \itemize{
#'   \item \code{H5Group} (CSC sparse) — calls \code{read_sparse_cols}
#'   \item \code{H5D} (dense) — reads via HDF5 hyperslab
#' }
#' When \code{col_idx} is \code{NULL} the existing \code{as.matrix} path is
#' used unchanged so the non-subsetting code path is unaffected.
#'
#' @param x An \code{H5Group} or \code{H5D} dataset representing a matrix
#'   slot (rows = features, columns = cells).
#' @param col_idx Integer vector of 1-based column (cell) indices, or
#'   \code{NULL} to read all columns.
#'
#' @return A matrix (dense or sparse as appropriate).
#'
#' @keywords internal
#'
read_assay_matrix <- function(x, col_idx = NULL) {
  if (is.null(col_idx)) {
    return(as.matrix(x = x))
  }
  if (inherits(x, 'H5Group')) {
    return(read_sparse_cols(grp = x, col_idx = col_idx))
  }
  # Dense H5D: use hyperslab to read only the requested columns
  x_dims <- x$dims
  return(as.matrix(x$read(args = list(seq_len(x_dims[1L]), col_idx))))
}

#' @return \code{dimnames}: returns a two-length list of character vectors for
#' row and column names. Row names should be in a column named \code{index}
#'
#' @aliases dimnames
#'
#' @rdname ReadH5
#' @method dimnames H5D
#' @export
#'
dimnames.H5D <- function(x) {
  if (!IsDType(x = x, dtype = 'H5T_COMPOUND')) {
    stop("'x' must be an HDF5 compound dataset", call. = FALSE)
  }
  colnames <- x$get_type()$get_cpd_labels()
  index <- match(x = 'index', table = colnames)
  if (is.na(x = index)) {
    rownames <- NULL
  } else {
    colnames <- colnames[-index]
    rownames <- unlist(
      x = x$read_low_level(mem_type = H5T_COMPOUND$new(
        labels = 'index',
        dtypes = x$get_type()$get_cpd_types()[[index]]
      )),
      use.names = FALSE
    )
  }
  return(list(rownames, colnames))
}
