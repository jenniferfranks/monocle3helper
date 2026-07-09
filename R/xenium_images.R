#' Build Xenium morphology image pixel data frames aligned to cell centroids
#'
#' Reads a single-channel morphology OME-TIFF pyramid level for each
#' Xenium sample (via reticulate + the Python \code{tifffile} package)
#' and converts it into a pixel-level data frame aligned to cell centroid
#' coordinates in a monocle3 cell_data_set.
#'
#' Xenium morphology images are typically tens of thousands of pixels per
#' side at full resolution, which is far too large to represent as a
#' pixel-per-row data.frame. Use \code{level > 0} to load a downsampled
#' pyramid level; this function assumes each pyramid level downsamples by
#' a factor of \code{2^level} relative to full resolution (override with
#' \code{downsample_factor} if a sample's pyramid uses a different scheme).
#'
#' @param cds A monocle3 cell_data_set with \code{x_centroid}/\code{y_centroid}
#'   and \code{sample} columns in colData (see \code{build_xenium_cds()}).
#' @param sample_table A data.frame produced by discover_xenium_sample_table().
#' @param channel Integer morphology channel index, matching the
#'   \code{morphology_focus_XXXX.ome.tif} file name (0 = DAPI by default).
#' @param level Pyramid level to load (0 = full resolution). Use a level
#'   high enough that the resulting image is tractable as a data.frame.
#' @param pixel_size_um Physical size, in microns, of one full-resolution
#'   (\code{level = 0}) pixel. Defaults to the standard Xenium value.
#' @param downsample_factor Optional override for the pyramid downsample
#'   factor at \code{level}. Defaults to \code{2^level}.
#' @param img_buffer Pixel buffer around cells used to flag image regions
#'   containing tissue.
#' @param conda_env Optional conda environment name containing
#'   \code{tifffile} (and \code{imagecodecs} for compressed OME-TIFFs).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{image}{Pixel-level data.frame with columns x, y, rgb.val,
#'       sample, sample_label, sample_group, contains_cells}
#'     \item{channel}{Channel index used}
#'     \item{level}{Pyramid level used}
#'     \item{pixel_size_um}{Full-resolution pixel size in microns}
#'     \item{downsample_factor}{Downsample factor applied at \code{level}}
#'     \item{effective_pixel_size_um}{Pixel size, in microns, at \code{level}}
#'   }
#'
#' @export
build_xenium_images <- function(
    cds,
    sample_table,
    channel           = 0,
    level             = 1,
    pixel_size_um     = 0.2125,
    downsample_factor = NULL,
    img_buffer        = 5,
    conda_env         = NULL
) {

  required_cols <- c("x_centroid", "y_centroid", "sample")

  if (!all(required_cols %in% colnames(SummarizedExperiment::colData(cds)))) {
    stop(
      "cds must contain cell centroid coordinates in colData: ",
      paste(required_cols, collapse = ", "),
      ". Did you run build_xenium_cds()?"
    )
  }

  if (!is.null(conda_env)) {
    reticulate::use_condaenv(conda_env, required = TRUE)
  }

  tifffile <- reticulate::import("tifffile")

  if (is.null(downsample_factor)) {
    downsample_factor <- 2^level
  }
  effective_pixel_size <- pixel_size_um * downsample_factor

  df_cells <- as.data.frame(SummarizedExperiment::colData(cds))
  df_cells$sample <- as.character(df_cells$sample)

  img_list <- list()

  for (i in seq_len(nrow(sample_table))) {

    sample_id <- as.character(sample_table$sample_id[i])
    df_samp   <- df_cells[df_cells$sample == sample_id, , drop = FALSE]

    if (nrow(df_samp) == 0) next

    morph_dir <- sample_table$morphology_dir[i]
    if (is.na(morph_dir) || !dir.exists(morph_dir)) next

    morph_file <- file.path(
      morph_dir,
      sprintf("morphology_focus_%04d.ome.tif", channel)
    )

    if (!file.exists(morph_file)) {
      warning("Morphology file not found, skipping sample ", sample_id,
              ": ", morph_file)
      next
    }

    message(
      "Reading Xenium morphology image: ", morph_file,
      " (channel ", channel, ", level ", level, ")"
    )

    img_py  <- tifffile$imread(morph_file, is_ome = FALSE, level = as.integer(level))
    img_mat <- as.matrix(img_py)

    n_pixels <- length(img_mat)
    if (n_pixels > 5e6) {
      warning(
        "Loaded image for sample ", sample_id, " has ", n_pixels,
        " pixels at level ", level, ". Consider a higher level for a ",
        "more tractable pixel data.frame."
      )
    }

    dims <- dim(img_mat)
    value <- as.vector(img_mat)
    value_range <- range(value, na.rm = TRUE)
    value_norm <- if (diff(value_range) > 0) {
      (value - value_range[1]) / diff(value_range)
    } else {
      rep(0, length(value))
    }

    df_img <- data.frame(
      y   = rep(seq_len(dims[1]), times = dims[2]),
      x   = rep(seq_len(dims[2]), each  = dims[1])
    )
    df_img$rgb.val <- grDevices::grey(value_norm)

    df_img$sample       <- sample_id
    df_img$sample_label <- sample_table$sample_label[i]
    df_img$sample_group <- sample_table$sample_group[i]

    x_rng <- range(df_samp$x_centroid / effective_pixel_size, na.rm = TRUE)
    y_rng <- range(df_samp$y_centroid / effective_pixel_size, na.rm = TRUE)

    df_img$contains_cells <- ifelse(
      df_img$x >= (x_rng[1] - img_buffer) &
        df_img$x <= (x_rng[2] + img_buffer) &
        df_img$y >= (y_rng[1] - img_buffer) &
        df_img$y <= (y_rng[2] + img_buffer),
      "yes", "no"
    )

    img_list[[length(img_list) + 1]] <- df_img
  }

  list(
    image                   = if (length(img_list) > 0) dplyr::bind_rows(img_list) else NULL,
    channel                 = channel,
    level                   = level,
    pixel_size_um           = pixel_size_um,
    downsample_factor       = downsample_factor,
    effective_pixel_size_um = effective_pixel_size
  )
}
